# Section 7.10 of FIMD
# https://stefvanbuuren.name/fimd/sec-mlguidelines.html

library(mice)
library(dplyr)
library(broom.mixed)
data("brandsma", package = "mice")
dat <- brandsma[, c("sch", "pup", "lpo",
                    "iqv", "ses", "ssi")]
md.pattern(dat)


# 7.10.1 ------------------------------------------------------------------


d <- brandsma[, c("sch", "lpo")]
pred <- make.predictorMatrix(d)
pred["lpo", "sch"] <- -2
imp <- mice(d, pred = pred, meth = "2l.pmm", m = 10, maxit = 1,
            print = FALSE, seed = 152)
library(lme4)
fit <- with(imp, lmer(lpo ~ (1 | sch), REML = FALSE))
library(broom.mixed)
summary(pool(fit))
library(mitml)
testEstimates(as.mitml.result(fit), extra.pars = TRUE)$extra.pars


# 7.10.2 ------------------------------------------------------------------


d <- brandsma[, c("sch", "lpo", "iqv")]
pred <- make.predictorMatrix(d)

# An entry of -2 in the predictor matrix signals the cluster variable, whereas
# an entry of 3 indicates that the cluster means of the covariates are added as
# a predictor to the imputation model.

pred["lpo", ] <- c(-2, 0, 3)
pred["iqv", ] <- c(-2, 3, 0)
pred

# lpo is imputed from iqv and the cluster means of iqv, while iqv is imputed
# from lpo and the cluster means of lpo. If the residuals are close to normal
# and the within-cluster error variances are similar, then 2l.pan is also a good
# choice.

imp <- mice(d, pred = pred, meth = "2l.pmm", seed = 919,
            m = 10, print = FALSE)
fit <- with(imp, lmer(lpo ~  iqv + (1 | sch), REML = FALSE))
summary(pool(fit))
testEstimates(as.mitml.result(fit), extra.pars = TRUE)$extra.pars


# 7.10.3 ------------------------------------------------------------------



library(dplyr)
res <- mice::complete(imp, "long")  |> 
  group_by(sch, .imp) |>
  mutate(iqm = mean(iqv)) |>
  group_by(.imp) |>
  do(model = lmer(lpo ~ iqv + iqm + (1 | sch),
                  REML = FALSE, data = .)) |>
  as.list() |> .[[-1]]
summary(pool(res))

# CF: updated code that doesn't use superceded `do()` function
res <- mice::complete(imp, "long")  |> 
  group_by(sch, .imp) |>
  mutate(iqm = mean(iqv)) |>
  ungroup() |> 
  nest_by(.imp) |>
  mutate(model = list(lmer(lpo ~ iqv + iqm + (1 | sch),
                  REML = FALSE, data = data))) |>
  as.list() |> 
  purrr::pluck("model")

summary(pool(res))

# 7.10.4

d <- brandsma[, c("sch", "lpo", "iqv", "den")]

# The missing values occur in lpo, iqv and den. den is a measured variable, so
# the value is identical for all members of the same cluster. If den is missing,
# it is missing for the entire cluster. Imputing a missing level-2 predictor is
# done by forming an imputation model at the cluster level.

meth <- make.method(d)
meth[c("lpo", "iqv", "den")] <- c("2l.pmm", "2l.pmm",
                                  "2lonly.pmm")
# 2lonly.pmm aggregates level-1 predictors and imputes the level-2 variable by
# predictive mean matching.

# Imputation methods for level-2 predictors should assign the same imputed value
# to all members within the same class

pred <- make.predictorMatrix(d)
pred["lpo", ] <- c(-2, 0, 3, 1)
pred["iqv", ] <- c(-2, 3, 0, 1)
pred["den", ] <- c(-2, 1, 1, 0)
imp <- mice(d, pred = pred, meth = meth, seed = 418,
            m = 10, print = FALSE)
# Fig 7.6
densityplot(imp, ~ lpo + den)

fit <- with(imp, lmer(lpo ~ 1 + iqv + as.factor(den)
                      + (1 | sch), REML = FALSE))
summary(pool(fit))
testEstimates(as.mitml.result(fit), extra.pars = TRUE)$extra.pars


# 7.10.5 ------------------------------------------------------------------



d <- brandsma[, c("sch", "lpo", "iqv", "sex", "den")]

# The new variables lpm, iqm and sxm will hold the cluster means of lpo, iqv and
# sex, respectively. Variables iqd and lpd will hold the values of iqv and lpo
# in deviations from their cluster means. Variables iqd.sex, lpd.sex and iqd.lpd
# are two-way interactions of level-1 variables scaled as deviations from the
# cluster means. Variables iqd.den, sex.den and lpd.den are cross-level
# interactions. Finally, iqm.den, sxm.den and lpm.den are interactions at
# level-2.

# The idea is that we impute lpo, iqv, sex and den, and update the other
# variables accordingly.
d <- data.frame(d, 
                lpm = NA, iqm = NA, sxm = NA,
                iqd = NA, lpd = NA,
                iqd.sex = NA, lpd.sex = NA, iqd.lpd = NA,
                iqd.den = NA, sex.den = NA, lpd.den = NA,
                iqm.den = NA, sxm.den = NA, lpm.den = NA)

# level-1 variables
meth <- make.method(d)
meth[c("lpo", "iqv", "sex")] <- "2l.pmm"

pred <- make.predictorMatrix(d)
pred[,] <- 0
pred[, "sch"] <- -2
codes <- c(3, 3, rep(1, 6))
pred["lpo", c("iqv", "sex", "iqd.sex", "sex.den", "iqd.den",
              "den", "iqm.den", "sxm.den")] <- codes
pred["iqv", c("lpo", "sex", "lpd.sex", "sex.den", "lpd.den",
              "den", "lpm.den", "sxm.den")] <- codes
pred["sex", c("lpo", "iqv", "iqd.lpd", "lpd.den", "iqd.den",
              "den", "iqm.den", "lpm.den")] <- codes

# level-2 variables
meth["den"] <- "2lonly.pmm"
pred["den", c("lpo", "iqv", "sex",
              "iqd.sex", "lpd.sex", "iqd.lpd")] <- 1

# The 2l.groupmean method from the miceadds package returns the cluster mean
# pertaining to each observation. Centering on the cluster means is widely
# practiced, but significantly alters the multilevel model. In the context of
# imputation, centering on the cluster means often enhances stability and
# robustness of models to generate imputations, especially if interactions are
# involved. When the complete-data model uses cluster centering, then the
# imputation model should also do so.

library(miceadds)
# derive group means
meth[c("iqm", "sxm", "lpm")] <- "2l.groupmean"
pred[c("iqm", "sxm", "lpm"), c("iqv", "sex", "lpo")] <- diag(3)

# derive deviations from cluster mean
meth["iqd"] <- "~ I(iqv - iqm)"
meth["lpd"] <- "~ I(lpo - lpm)"

# The next block of code specifies the interaction effects, by means of passive
# imputation.

# derive interactions
meth["iqd.sex"] <- "~ I(iqd * sex)"
meth["lpd.sex"] <- "~ I(lpd * sex)"
meth["iqd.lpd"] <- "~ I(iqd * lpd)"
meth["iqd.den"] <- "~ I(iqd * den)"
meth["sex.den"] <- "~ I(sex * den)"
meth["lpd.den"] <- "~ I(lpd * den)"
meth["iqm.den"] <- "~ I(iqm * den)"
meth["sxm.den"] <- "~ I(sxm * den)"
meth["lpm.den"] <- "~ I(lpm * den)"

# The visit sequence specified below updates the relevant derived variables
# after any of the measured variables is imputed, so that interactions are
# always in sync.

visit <- c("lpo", "lpm", "lpd",
           "lpd.sex", "iqd.lpd", "lpd.den", "lpm.den",
           "iqv", "iqm", "iqd",
           "iqd.sex", "iqd.lpd", "iqd.den", "iqm.den",
           "sex", "sxm",
           "iqd.sex", "lpd.sex", "sex.den", "sxm.den",
           "den", "iqd.den", "sex.den", "lpd.den",
           "iqm.den", "sxm.den", "lpm.den")

imp <- mice(d, pred = pred, meth = meth, seed = 188,
            visit = visit, m = 10, print = FALSE,
            allow.na = TRUE)

long <- mice::complete(imp, "long", include = TRUE)
long$den <- as.factor(long$den)
imp2 <- as.mids(long)

# three types of multiplicative interactions among the predictors into the
# model:
# a level-1 interaction: iqv*sex (both change within school)
# a cross-level interaction: sex*den (sex changes within school, den does not)
# level-2 interaction: iqm*den (both do not change within school)

fit <- with(imp2, lmer(lpo ~ 1 + iqv*sex + iqm*den + sex*den
                       + (1 | sch), REML = FALSE))
summary(pool(fit))


# 7.10.6 ------------------------------------------------------------------



d <- brandsma[, c("sch", "lpo", "iqv")]

# lpo needs to be centered around the grand mean in order to reduce the large
# number of warnings about unstable estimates.

d$lpo <- as.vector(scale(d$lpo, scale = FALSE)) # center
pred <- make.predictorMatrix(d)

# The entry of 4 at cell (lpo, iqv) in the predictor matrix adds three variables
# to the imputation model for lpo: the value of iqv, the cluster means of iqv
# and the random slopes of iqv. Conversely, imputing iqv adds the three
# covariates: the values of lpo, the cluster means of lpo and the random slopes
# of lpo.

pred["lpo", ] <- c(-2, 0, 4)
pred["iqv", ] <- c(-2, 4, 0)
pred

imp <- mice(d, pred = pred, meth = "2l.pmm", seed = 441,
            m = 10, print = FALSE, maxit = 20)

# Grand-mean centering implies a little extra work because we must
# back-transform the data if we want the values in the original scale. What
# remains is that rescaling improves speed and stability, so for the purpose of
# imputation I recommend to scale level-1 variables in deviations from their
# means.

# The following code block unfolds the mids object, adds the IQ cluster means,
# restores the rescaling of lpo, and estimates and combines the parameters of
# the random slopes model.

imp2 <- mice::complete(imp, "long", include = TRUE) |>
  group_by(sch) |>
  mutate(iqm = mean(iqv, na.rm = TRUE),
         lpo = lpo + mean(brandsma$lpo, na.rm = TRUE)) |>
  as.mids()
fit <- with(imp2, lmer(lpo ~  iqv + iqm + (1 + iqv | sch),
                       REML = FALSE))
summary(pool(fit))
testEstimates(as.mitml.result(fit), extra.pars = TRUE)$extra.pars


# 7.10.7 ------------------------------------------------------------------

# CF: this section is incomplete

d <- brandsma[, c("sch", "lpo", "iqv", "ses")]
d$lpo <- as.vector(scale(d$lpo, scale = FALSE))
d <- data.frame(d,
                iqv.ses = NA, ses.lpo = NA, iqv.lpo = NA,
                lpm = NA, iqm = NA, sem = NA,
                iqv.iqm = NA, ses.sem = NA, lpo.lpm = NA,
                iqv.sem = NA, iqv.lpm = NA,
                ses.iqm = NA, ses.lpm = NA,
                lpo.iqm = NA, lpo.sem = NA,
                iqm.sem = NA, lpm.sem = NA, iqm.lpm = NA)

meth <- make.method(d)
meth[c("iqm", "sem", "lpm")] <- "2l.groupmean"
meth["iqv.ses"] <- "~ I(iqv * ses)"
meth["ses.lpo"] <- "~ I(ses * lpo)"
meth["iqv.lpo"] <- "~ I(iqv * lpo)"

meth["iqv.iqm"] <- "~ I(iqv * lpo)"
meth["ses.sem"] <- "~ I(ses * sem)"
meth["lpo.lpm"] <- "~ I(lpo * lpm)"
meth["iqv.sem"] <- "~ I(iqv * sem)"
meth["iqv.lpm"] <- "~ I(iqv * lpm)"
meth["ses.iqm"] <- "~ I(ses * iqm)"
meth["ses.lpm"] <- "~ I(ses * lpm)"
meth["lpo.iqm"] <- "~ I(lpo * iqm)"
meth["lpo.sem"] <- "~ I(lpo * sem)"
meth["iqm.sem"] <- "~ I(iqm * sem)"
meth["lpm.sem"] <- "~ I(lpm * sem)"
meth["iqm.lpm"] <- "~ I(iqm *lpm)"


pred <- make.predictorMatrix(d)
pred[c("iqm", "sem", "lpm"), ] <- 0
pred["iqm", c("sch", "iqv")] <- c(-2, 1)
pred["sem", c("sch", "ses")] <- c(-2, 1)
pred["lpm", c("sch", "lpo")] <- c(-2, 1)

visit <- c("lpo", "iqv.lpo", "ses.lpo",
           "lpm", "lpo.lpm", "iqv.lpm", "ses.lpm",
           "lpo.iqm", "lpo.sem", "iqm.lpm", "lpm.sem",
           "iqv", "iqv.ses", "iqv.lpo",
           "iqm", "iqv.iqm", "iqv.sem", "iqv.lpm",
           "ses.iqm", "lpo.iqm", "iqm.sem", "iqm.lpm",
           "ses", "iqv.ses", "ses.lpo",
           "sem", "ses.sem", "iqv.sem", "ses.iqm",
           "ses.lpm", "lpo.sem", "iqm.sem", "lpm.sem")

imp <- mice(d, pred = pred, meth = meth, seed = 211,
            visit = visit, m = 10, print = FALSE, maxit = 10,
            allow.na = TRUE)
fit <- with(imp, lmer(lpo ~ iqv * ses + iqm * sem +
                        iqv * iqm + iqv * sem +
                        ses * iqm + ses * sem + (1 + ses + iqv | sch),
                      REML = FALSE))
summary(pool(fit))
testEstimates(as.mitml.result(fit), extra.pars = TRUE)$extra.pars

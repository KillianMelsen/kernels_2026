library(asreml)
asreml.options(workspace = "5000mb", pworkspace = "5000mb", maxit = 50)
set.seed(1997)
d <- readRDS("DROPS/data/ydata.SE.rds")
K <- readRDS("DROPS/data/K.SE.rds")[levels(d$Gen), levels(d$Gen)]

# Fitting GY ~ -1 + Man + Env + Man:Env is identical to fitting GY ~ -1 + ManEnv!
fit1 <- asreml(fixed = GY ~ -1 + Man + Env + Man:Env,
               random = ~ vm(Gen, K) + idv(Man):vm(Gen, K) + idv(Env):vm(Gen, K) + idv(ManEnv):vm(Gen, K),
               residual = ~ units,
               data = d,
               weights = weight,
               family = asr_gaussian(dispersion = 1.0))

summary(fit1)$varcomp

# Load in results from variance partitioning:
r <- readRDS("DROPS_LOF/results_DROPS_LOF_SE.rds")
aggregate(r, cbind(Vge, Vlof, Ve) ~ Model, FUN = mean)

# numbers for the model with all main effects and (lower-order) interactions:
c(summary(fit1)$varcomp[1, 1], #.......Variance captured by the random genetic main effect, this is like fitting a GxE term with all correlations equal to 1.
  sum(summary(fit1)$varcomp[2, 1], #...Variance captured by the random management x genotype interaction, first component of the lack of fit variance.
      summary(fit1)$varcomp[3, 1], #...Variance captured by the random environment x genotype interaction, second component of the lack of fit variance.
      summary(fit1)$varcomp[4, 1]), #..Variance captured by the random management x environment x genotype interaction, third component of the lack of fit variance.
  mean(1 / d$weight)) #................Average residual variance. Weights are the diagonal elements of the precision matrices of BLUEs from the single-trial models fitted using LMMsolver.

# Compare that to the numbers from the lack of fit model (see ADD model in /DROPS_LOF/DROPS_LOF.r) which has a fixed ManEnv effect, a random genetic main effect, and a random idh(ManEnv):vm(Gen, K) three-way interaction.
aggregate(r, cbind(Vge, Vlof, Ve) ~ Model, FUN = mean)[1, 2:4]

# Let's make things more complicated:
fit2 <- asreml(fixed = GY ~ -1 + Man + Env + Man:Env,
               random = ~ vm(Gen, K) + idh(Man):vm(Gen, K) + idh(Env):vm(Gen, K) + idh(ManEnv):vm(Gen, K),
               residual = ~ units,
               data = d,
               weights = weight,
               family = asr_gaussian(dispersion = 1.0))

summary(fit2)$varcomp

# Numbers according to this more complicated model:
c(summary(fit2)$varcomp[1, 1], #.................Variance captured by the random genetic main effect, this is like fitting a GxE term with all correlations equal to 1.
  sum(mean(summary(fit2)$varcomp[2:3, 1]), #.....Average variance captured by the random management x genotype interaction, first component of the lack of fit variance.
      mean(summary(fit2)$varcomp[4:17, 1]), #....Average variance captured by the random environment x genotype interaction, second component of the lack of fit variance.
      mean(summary(fit2)$varcomp[18:45, 1])), #..Average variance captured by the random management x environment x genotype interaction, third component of the lack of fit variance.
  mean(1 / d$weight)) #..........................Average residual variance. Weights are the diagonal elements of the precision matrices of BLUEs from the single-trial models fitted using LMMsolver.




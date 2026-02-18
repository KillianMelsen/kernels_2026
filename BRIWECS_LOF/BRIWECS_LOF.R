library(asreml)
library(ggplot2)
asreml.options(workspace = "5000mb", pworkspace = "5000mb", maxit = 50)
set.seed(1997)
envs <- levels(readRDS("Briwecs/data/datalist.rds")$ydata$Env)

# Loading kinship, env. correlation matrix and full data:
d.full <- readRDS("Briwecs/data/datalist.rds")$ydata
d.full <- droplevels(d.full[d.full$Env %in% envs, c("Env", "Man", "Gen", "GY")])
d.full$EnvGen <- as.factor(paste(d.full$Env, d.full$Gen, sep = ":"))
d.full$ManGen <- as.factor(paste(d.full$Man, d.full$Gen, sep = ":"))
d.full$ManEnv <- as.factor(paste(d.full$Man, d.full$Env, sep = ":"))
d.full$ManEnvGen <- as.factor(paste(d.full$Man, d.full$Env, d.full$Gen, sep = ":"))
d.full$ManEnv2 <- d.full$ManEnv

# Kinship and environmental matrices:
K <- readRDS("Briwecs/data/K.rds")[levels(d.full$Gen), levels(d.full$Gen)]
EC <- readRDS("Briwecs/data/datalist.rds")$EC[levels(d.full$Env), levels(d.full$Env)]
ED <- readRDS("Briwecs/data/datalist.rds")$ED[levels(d.full$Env), levels(d.full$Env)]

# Results structure ====
results <- expand.grid(Model = c("FA-1", "FA-2", "FA-3", "SV-LK", "SV-GK", "MV-LK", "MV-GK"),
                       Management = levels(d.full$Man),
                       Environment = levels(d.full$Env),
                       Vge = NA,
                       Vlof = NA,
                       Ve = NA)

# FA1 ====
mod.FA1 <- asreml(GY ~ -1 + ManEnv,
                  random = ~ fa(ManEnv, 1):vm(Gen, K),
                  residual = ~ units,
                  data = d.full,
                  trace = TRUE)

PSI <- summary(mod.FA1)$varcomp[1:28, "component"]
names(PSI) <- gsub(".*!(.*)!var", "\\1", rownames(summary(mod.FA1)$varcomp)[1:28])
L <- matrix(0, 28, 1)
rownames(L) <- gsub(".*!(.*)!var", "\\1", rownames(summary(mod.FA1)$varcomp)[1:28])
L[, 1] <- summary(mod.FA1)$varcomp[29:56, "component"]
LLvars <- diag(L %*% t(L))
for (i in 1:28) {
  man <- gsub("(..):(.*)", "\\1", names(PSI)[i])
  env <- gsub("(..):(.*)", "\\2", names(PSI)[i])
  results[results$Model == "FA-1" & results$Management == man & results$Environment == env, "Vge"] <- LLvars[i]
  results[results$Model == "FA-1" & results$Management == man & results$Environment == env, "Vlof"] <- PSI[i]
}
results[results$Model == "FA-1", "Ve"] <- summary(mod.FA1)$varcomp[nrow(summary(mod.FA1)$varcomp), "component"]

# FA2 ====
mod.FA2 <- asreml(GY ~ -1 + ManEnv,
                  random = ~ fa(ManEnv, 2):vm(Gen, K),
                  residual = ~ units,
                  data = d.full,
                  trace = TRUE)

PSI <- summary(mod.FA2)$varcomp[1:28, "component"]
names(PSI) <- gsub(".*!(.*)!var", "\\1", rownames(summary(mod.FA2)$varcomp)[1:28])
L <- matrix(0, 28, 2)
rownames(L) <- gsub(".*!(.*)!var", "\\1", rownames(summary(mod.FA2)$varcomp)[1:28])
L[, 1] <- summary(mod.FA2)$varcomp[29:56, "component"]
L[, 2] <- summary(mod.FA2)$varcomp[57:84, "component"]
LLvars <- diag(L %*% t(L))
for (i in 1:28) {
  man <- gsub("(..):(.*)", "\\1", names(PSI)[i])
  env <- gsub("(..):(.*)", "\\2", names(PSI)[i])
  results[results$Model == "FA-2" & results$Management == man & results$Environment == env, "Vge"] <- LLvars[i]
  results[results$Model == "FA-2" & results$Management == man & results$Environment == env, "Vlof"] <- PSI[i]
}
results[results$Model == "FA-2", "Ve"] <- summary(mod.FA2)$varcomp[nrow(summary(mod.FA2)$varcomp), "component"]

# FA3 ====
mod.FA3 <- asreml(GY ~ -1 + ManEnv,
                  random = ~ fa(ManEnv, 3):vm(Gen, K),
                  residual = ~ units,
                  data = d.full,
                  trace = TRUE)

PSI <- summary(mod.FA3)$varcomp[1:28, "component"]
names(PSI) <- gsub(".*!(.*)!var", "\\1", rownames(summary(mod.FA3)$varcomp)[1:28])
L <- matrix(0, 28, 3)
rownames(L) <- gsub(".*!(.*)!var", "\\1", rownames(summary(mod.FA3)$varcomp)[1:28])
L[, 1] <- summary(mod.FA3)$varcomp[29:56, "component"]
L[, 2] <- summary(mod.FA3)$varcomp[57:84, "component"]
L[, 3] <- summary(mod.FA3)$varcomp[85:112, "component"]
LLvars <- diag(L %*% t(L))
for (i in 1:28) {
  man <- gsub("(..):(.*)", "\\1", names(PSI)[i])
  env <- gsub("(..):(.*)", "\\2", names(PSI)[i])
  results[results$Model == "FA-3" & results$Management == man & results$Environment == env, "Vge"] <- LLvars[i]
  results[results$Model == "FA-3" & results$Management == man & results$Environment == env, "Vlof"] <- PSI[i]
}
results[results$Model == "FA-3", "Ve"] <- summary(mod.FA3)$varcomp[nrow(summary(mod.FA3)$varcomp), "component"]

# Single-var relmat model version 2 ====
vf <- function(order, kappa) {
  # kappa[1] = variance for all Env within trait 1
  # kappa[2] = variance for all Env within trait 2
  # kappa[3] = correlation between M1 and M2
  # The correlation matrix of the Man levels:
  Rm <- matrix(1, 2, 2)
  Rm[1, 2] <- Rm[2, 1] <- kappa[3]
  
  # The full covariance matrix:
  S <- outer(sqrt(kappa[1:ncol(Rm)]), sqrt(kappa[1:ncol(Rm)]))
  V <- kronecker(S * Rm, EC)
  
  # Derivatives wrt kappa[1] and kappa[2] (variances)
  varderivs <- vector("list", ncol(Rm))
  for (dk in 1:ncol(Rm)) {
    # Indicator matrix of where kappa[dk] is present:
    I <- matrix(0, nrow(Rm), ncol(Rm))
    I[dk,] <- I[, dk] <- 1
    I <- kronecker(I, matrix(1, nrow(EC), ncol(EC)))
    tmp <- sqrt(kappa[1:ncol(Rm)])
    tmp[dk] <- 1 / tmp[dk]
    tmp <- outer(tmp, tmp)
    tmp[dk, dk] <- 1
    deriv <- 0.5 * I * kronecker(tmp * Rm, EC)
    deriv[((dk - 1) * nrow(EC) + 1):(dk * nrow(EC)), ((dk - 1) * nrow(EC) + 1):(dk * nrow(EC))] <-
      deriv[((dk - 1) * nrow(EC) + 1):(dk * nrow(EC)), ((dk - 1) * nrow(EC) + 1):(dk * nrow(EC))] * 2
    varderivs[[dk]] <- deriv
  }
  
  # Derivative wrt kappa[3]
  # Indicator matrix of where kappa[3] is present:
  I <- matrix(1, nrow(Rm), ncol(Rm))
  diag(I) = 0
  dkcorr <- kronecker(S * I, EC)
  
  return(c(list(V), varderivs, list(dkcorr)))
}

init <- c(0.1, 0.1, 0.1)
type <- c("V", "V", "R")
con <- c("P", "P", "U")
mod.svar.EC <- asreml(fixed = GY ~ -1 + ManEnv,
                      random = ~ own(ManEnv, "vf", init, type, con):vm(Gen, K) + diag(ManEnv2):vm(Gen, K),
                      residual = ~ units,
                      data = d.full,
                      trace = TRUE)

results[results$Model == "SV-LK" & results$Management == "HN", "Vge"] <- summary(mod.svar.EC)$varcomp[1, "component"]
results[results$Model == "SV-LK" & results$Management == "LN", "Vge"] <- summary(mod.svar.EC)$varcomp[2, "component"]
results[results$Model == "SV-LK", "Ve"] <- summary(mod.svar.EC)$varcomp[nrow(summary(mod.svar.EC)$varcomp), "component"]
for (i in 4:(nrow(summary(mod.svar.EC)$varcomp) - 1)) {
  management <- gsub(".*_([HL]N):(.*)", "\\1", rownames(summary(mod.svar.EC)$varcomp)[i])
  environment <- gsub(".*_([HL]N):(.*)", "\\2", rownames(summary(mod.svar.EC)$varcomp)[i])
  results[results$Model == "SV-LK" & results$Management == management & results$Environment == environment, "Vlof"] <- summary(mod.svar.EC)$varcomp$component[i]
}

# Single-var Gaussian kernel model 2 ====
vf <- function(order, kappa) {
  # kappa[1] = variance of all Env within trait/M 1
  # kappa[2] = variance of all Env within trait/M 2
  # kappa[3] = correlation between M1 and M2
  # kappa[4] = bandwidth parameter of the Gaussian kernel for Env
  # The correlation matrix of the Man levels:
  Rm <- matrix(1, 2, 2)
  Rm[1, 2] <- Rm[2, 1] <- kappa[3]
  
  # The full covariance matrix:
  S <- outer(sqrt(kappa[1:ncol(Rm)]), sqrt(kappa[1:ncol(Rm)]))
  V <- kronecker(S * Rm, exp(-kappa[4] * ED))
  
  # Derivatives wrt kappa[1] and kappa[2] (variances)
  varderivs <- vector("list", ncol(Rm))
  for (dk in 1:ncol(Rm)) {
    # Indicator matrix of where kappa[dk] is present:
    I <- matrix(0, nrow(Rm), ncol(Rm))
    I[dk,] <- I[, dk] <- 1
    I <- kronecker(I, matrix(1, nrow(ED), ncol(ED)))
    tmp <- sqrt(kappa[1:ncol(Rm)])
    tmp[dk] <- 1 / tmp[dk]
    tmp <- outer(tmp, tmp)
    tmp[dk, dk] <- 1
    deriv <- 0.5 * I * kronecker(tmp * Rm, exp(-kappa[4] * ED))
    deriv[((dk - 1) * nrow(ED) + 1):(dk * nrow(ED)), ((dk - 1) * nrow(ED) + 1):(dk * nrow(ED))] <-
      deriv[((dk - 1) * nrow(ED) + 1):(dk * nrow(ED)), ((dk - 1) * nrow(ED) + 1):(dk * nrow(ED))] * 2
    varderivs[[dk]] <- deriv
  }
  
  # Derivative wrt kappa[3]
  # Indicator matrix of where kappa[2] is present in Rm:
  I <- matrix(1, nrow(Rm), ncol(Rm))
  diag(I) = 0
  dkcorr <- kronecker(S * I, exp(-kappa[4] * ED))
  
  # Derivative wrt kappa[4]
  # Indicator matrix of where kappa[4] is present in e^(-h*ED)
  IRm <- matrix(1, nrow(ED), ncol(ED))
  IED <- matrix(1, nrow(Rm), ncol(Rm))
  dkbw <- (kronecker(S * Rm, IRm) * kronecker(IED, -ED) * kronecker(IED, exp(-kappa[4] * ED)))
  
  # cat(kappa, "\n\n")
  
  return(c(list(V), varderivs, list(dkcorr, dkbw)))
}

init <- c(0.1, 0.1, 0.1, 0.1)
type <- c("V", "V", "R", "V")
con <- c("P", "P", "U", "P")
mod.svar.GK <- asreml(fixed = GY ~ -1 + ManEnv,
                      random = ~ own(ManEnv, "vf", init, type, con):vm(Gen, K) + diag(ManEnv2):vm(Gen, K),
                      residual = ~ units,
                      data = d.full,
                      trace = TRUE)

results[results$Model == "SV-GK" & results$Management == "HN", "Vge"] <- summary(mod.svar.GK)$varcomp[1, "component"]
results[results$Model == "SV-GK" & results$Management == "LN", "Vge"] <- summary(mod.svar.GK)$varcomp[2, "component"]
results[results$Model == "SV-GK", "Ve"] <- summary(mod.svar.GK)$varcomp[nrow(summary(mod.svar.GK)$varcomp), "component"]
for (i in 5:(nrow(summary(mod.svar.GK)$varcomp) - 1)) {
  management <- gsub(".*_([HL]N):(.*)", "\\1", rownames(summary(mod.svar.GK)$varcomp)[i])
  environment <- gsub(".*_([HL]N):(.*)", "\\2", rownames(summary(mod.svar.GK)$varcomp)[i])
  results[results$Model == "SV-GK" & results$Management == management & results$Environment == environment, "Vlof"] <- summary(mod.svar.GK)$varcomp$component[i]
}

# Multi-var relmat model 2 ====
vf <- function(order, kappa) {
  # kappa[1] = variance of M1:E1
  # kappa[2] = variance of M1:E2
  #   ...
  # kappa[p * q] = variance of Mp:Eq
  # kappa[p * q + 1] = correlation between M1 and M2
  # Number of managements:
  n.mans <- order / nrow(EC)
  
  # The correlation matrix of the Man levels (specify manually!):
  Rm <- matrix(1, 2, 2)
  Rm[1, 2] <- Rm[2, 1] <- kappa[order + 1]
  
  # The full covariance matrix:
  S <- outer(sqrt(kappa[1:order]), sqrt(kappa[1:order]))
  V <- S * kronecker(Rm, EC)
  
  # Derivative wrt kappa[p * q + 1]
  # Indicator matrix of where kappa[p*q+1] is present:
  I <- matrix(1, nrow(Rm), ncol(Rm))
  diag(I) = 0
  I <- kronecker(I, matrix(1, nrow(EC), ncol(EC)))
  dkrm <- (S * I) * kronecker(matrix(1, nrow(Rm), ncol(Rm)), EC)
  
  # Derivatives wrt all variances
  varderivs <- vector("list", order)
  for (dk in 1:order) {
    # Indicator matrix of where kappa[dk] is present:
    I <- matrix(0, order, order)
    I[dk,] <- I[, dk] <- 1
    tmp <- sqrt(kappa[1:order])
    tmp[dk] <- 1 / tmp[dk]
    tmp <- outer(tmp, tmp)
    tmp[dk, dk] <- 1
    deriv <- 0.5 * I * tmp * kronecker(Rm, EC)
    deriv[dk, dk] <- 1
    varderivs[[dk]] <- deriv
  }
  # cat(kappa, "\n\n")
  return(c(list(V), varderivs, list(dkrm)))
}

init <- c(rep(0.1, nrow(EC) * length(levels(d.full$Man))), 0.1)
type <- c(rep("V", nrow(EC) * length(levels(d.full$Man))), "R")
con <- c(rep("P", nrow(EC) * length(levels(d.full$Man))), "U")
mod.mvar.EC <- asreml(fixed = GY ~ -1 + ManEnv,
                      random = ~ own(ManEnv, "vf", init, type):vm(Gen, K) + diag(ManEnv2):vm(Gen, K),
                      residual = ~ units,
                      data = d.full,
                      trace = TRUE)

genvars <- summary(mod.mvar.EC)$varcomp[1:28, "component"]
names(genvars) <- levels(d.full$ManEnv)
for (i in 1:28) {
  man <- gsub("(..):(.*)", "\\1", names(genvars)[i])
  env <- gsub("(..):(.*)", "\\2", names(genvars)[i])
  results[results$Model == "MV-LK" & results$Management == man & results$Environment == env, "Vge"] <- genvars[i]
}
results[results$Model == "MV-LK", "Ve"] <- summary(mod.mvar.EC)$varcomp[nrow(summary(mod.mvar.EC)$varcomp), "component"]
for (i in 30:(nrow(summary(mod.mvar.EC)$varcomp) - 1)) {
  management <- gsub(".*_([HL]N):(.*)", "\\1", rownames(summary(mod.mvar.EC)$varcomp)[i])
  environment <- gsub(".*_([HL]N):(.*)", "\\2", rownames(summary(mod.mvar.EC)$varcomp)[i])
  results[results$Model == "MV-LK" & results$Management == management & results$Environment == environment, "Vlof"] <- summary(mod.mvar.EC)$varcomp$component[i]
}

# Multi-var Gaussian kernel model 2 ====
vf <- function(order, kappa) {
  # kappa[1] = variance of M1:E1
  # kappa[2] = variance of M1:E2
  #   ...
  # kappa[p * q] = variance of Mp:Eq
  # kappa[p * q + 1] = correlation between M1 and M2
  # kappa[p * q + 2] = bandwidth parameter of the Gaussian kernel for Env
  # Number of managements:
  n.mans <- order / nrow(ED)
  
  # The correlation matrix of the Man levels (specify manually!):
  Rm <- matrix(1, 2, 2)
  Rm[1, 2] <- Rm[2, 1] <- kappa[order + 1]
  
  # The full covariance matrix:
  S <- outer(sqrt(kappa[1:order]), sqrt(kappa[1:order]))
  V <- S * kronecker(Rm, exp(-kappa[order + 2] * ED))
  
  # Derivative wrt kappa[p * q + 1]
  # Indicator matrix of where kappa[p * q + 1] is present:
  I <- matrix(1, nrow(Rm), ncol(Rm))
  diag(I) = 0
  I <- kronecker(I, matrix(1, nrow(ED), ncol(ED)))
  dkrm <- (S * I) * kronecker(matrix(1, nrow(Rm), ncol(Rm)), exp(-kappa[order + 2] * ED))
  
  # Derivative wrt kappa[p * q + 2]
  dkh <- S * kronecker(Rm, -ED * exp(-kappa[order + 2] * ED))
  
  # Derivatives wrt all variances
  varderivs <- vector("list", order)
  for (dk in 1:order) {
    # Indicator matrix of where kappa[dk] is present:
    I <- matrix(0, order, order)
    I[dk,] <- I[, dk] <- 1
    tmp <- sqrt(kappa[1:order])
    tmp[dk] <- 1 / tmp[dk]
    tmp <- outer(tmp, tmp)
    tmp[dk, dk] <- 1
    deriv <- 0.5 * I * tmp * kronecker(Rm, exp(-kappa[order + 2] * ED))
    deriv[dk, dk] <- 1
    varderivs[[dk]] <- deriv
  }
  # cat(kappa, "\n\n")
  return(c(list(V), varderivs, list(dkrm, dkh)))
}

init <- c(rep(0.1, nrow(ED) * length(levels(d.full$Man))), 0.1, 0.1)
type <- c(rep("V", nrow(ED) * length(levels(d.full$Man))), "R", "V")
con <- c(rep("P", nrow(ED) * length(levels(d.full$Man))), "U", "P")
mod.mvar.GK <- asreml(fixed = GY ~ -1 + ManEnv,
                      random = ~ own(ManEnv, "vf", init, type, con):vm(Gen, K) + diag(ManEnv2):vm(Gen, K),
                      residual = ~ units,
                      data = d.full,
                      trace = TRUE)

genvars <- summary(mod.mvar.GK)$varcomp[1:28, "component"]
names(genvars) <- levels(d.full$ManEnv)
for (i in 1:28) {
  man <- gsub("(..):(.*)", "\\1", names(genvars)[i])
  env <- gsub("(..):(.*)", "\\2", names(genvars)[i])
  results[results$Model == "MV-GK" & results$Management == man & results$Environment == env, "Vge"] <- genvars[i]
}
results[results$Model == "MV-GK", "Ve"] <- summary(mod.mvar.GK)$varcomp[nrow(summary(mod.mvar.GK)$varcomp), "component"]
for (i in 31:(nrow(summary(mod.mvar.GK)$varcomp) - 1)) {
  management <- gsub(".*_([HL]N):(.*)", "\\1", rownames(summary(mod.mvar.GK)$varcomp)[i])
  environment <- gsub(".*_([HL]N):(.*)", "\\2", rownames(summary(mod.mvar.GK)$varcomp)[i])
  results[results$Model == "MV-GK" & results$Management == management & results$Environment == environment, "Vlof"] <- summary(mod.mvar.GK)$varcomp$component[i]
}
  
saveRDS(results, "BRIWECS_LOF/results_BRIWECS_LOF.rds")
results <- readRDS("BRIWECS_LOF/results_BRIWECS_LOF.rds")
library(ggplot2)

# Plotting
colnames(results) <- c("Model", "Management", "Environment", "Covariables", "LOF", "Residual")
results2 <- as.data.frame(tidyr::pivot_longer(results, 4:6, names_to = "Component", values_to = "Variance"))
results3 <- as.data.frame(tidyr::pivot_longer(aggregate(results, cbind(Covariables, LOF, Residual) ~ Model + Management, FUN = mean), 3:5, names_to = "Component", values_to = "Variance"))
results2$Component <- factor(results2$Component, levels = c("Covariables", "LOF", "Residual"), labels = c("Covariables", "Lack of fit", "Residual"))
results3$Component <- factor(results3$Component, levels = c("Covariables", "LOF", "Residual"), labels = c("Covariables", "Lack of fit", "Residual"))
results2 <- droplevels(results2[results2$Model %in% c("FA-1", "FA-2", "FA-3", "SV-LK", "SV-GK", "MV-LK", "MV-GK"),])
results3 <- droplevels(results3[results3$Model %in% c("FA-1", "FA-2", "FA-3", "SV-LK", "SV-GK", "MV-LK", "MV-GK"),])

levels(results2$Management) <- c("High Nitrogen", "Low Nitrogen")
levels(results3$Management) <- c("High Nitrogen", "Low Nitrogen")

ggplot(results2, aes(fill = Component, y = Variance, x = Model)) +
  facet_grid(rows = vars(Environment), cols = vars(Management)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = c("#fcdd06", "#db161f", "#0e44af")) +
  scale_y_continuous(labels = scales::percent, breaks = c(0.0, 0.5, 1.0)) +
  ylab("Percentage of total variance") +
  theme_classic(base_size = 18) + theme(legend.position = "bottom")
ggsave(filename = "plots/BRIWECS_LOF.png", dpi = 300, width = 32, height = 40, units = "cm")

ggplot(results3, aes(fill = Component, y = Variance, x = Model)) +
  facet_grid(cols = vars(Management)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = c("#fcdd06", "#db161f", "#0e44af")) +
  scale_y_continuous(labels = scales::percent, breaks = c(0.0, 0.5, 1.0)) +
  ylab("Percentage of total\nvariance") +
  theme_classic(base_size = 18) + theme(legend.position = "bottom")
ggsave(filename = "plots/BRIWECS_LOF_Averaged.png", dpi = 300, width = 32, height = 10, units = "cm")

results$Total <- results$Covariables + results$LOF + results$Residual
results$GxExM_percentage <- results$Covariables / results$Total
results$LOF_percentage <- results$LOF / results$Total
results$Residual_percentage <- results$Residual / results$Total

# Output ====
## Figure 4 ====
ggplot(results3, aes(fill = Component, y = Variance, x = Model)) +
  facet_grid(cols = vars(Management)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#fcdd06", "#db161f", "#0e44af")) +
  # scale_y_continuous(labels = scales::percent, breaks = c(0.0, 0.5, 1.0)) +
  ylab("Variance") +
  theme_classic(base_size = 18) + theme(legend.position = "bottom") +
  ylim(c(0, 115))
ggsave(filename = "plots/BRIWECS_LOF_Averaged_numeric.png", dpi = 300, width = 32, height = 15, units = "cm")

ggplot(droplevels(results2[results2$Environment %in% levels(results2$Environment)[1:7],]), aes(fill = Component, y = Variance, x = Model)) +
  facet_grid(cols = vars(Management), rows = vars(Environment), scales = "free_y") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#fcdd06", "#db161f", "#0e44af")) +
  # scale_y_continuous(labels = scales::percent, breaks = c(0.0, 0.5, 1.0)) +
  ylab("Variance") +
  theme_classic(base_size = 18) + theme(legend.position = "bottom")
ggsave(filename = "plots/BRIWECS_LOF_perEnv_numeric_A.png", dpi = 300, width = 32, height = 48, units = "cm")

ggplot(droplevels(results2[results2$Environment %in% levels(results2$Environment)[8:14],]), aes(fill = Component, y = Variance, x = Model)) +
  facet_grid(cols = vars(Management), rows = vars(Environment), scales = "free_y") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#fcdd06", "#db161f", "#0e44af")) +
  # scale_y_continuous(labels = scales::percent, breaks = c(0.0, 0.5, 1.0)) +
  ylab("Variance") +
  theme_classic(base_size = 18) + theme(legend.position = "bottom")
ggsave(filename = "plots/BRIWECS_LOF_perEnv_numeric_B.png", dpi = 300, width = 32, height = 48, units = "cm")

## Text of section 3.2.2 ====
# Total genetic variance per management:
tmp <- aggregate(results, cbind(Covariables, LOF, Residual, Total) ~ Management, FUN = mean)
tmp$GenTotal <- round(tmp$Covariables + tmp$LOF, 2)
tmp[, c("Management", "GenTotal")]

# Percentages of genetic variance explained by the environmental covariables per model:
tmp <- aggregate(results, cbind(Covariables, LOF, Residual, Total) ~ Model, FUN = mean)
tmp$PercCov <- round(tmp$Covariables / (tmp$Covariables + tmp$LOF), 2)
tmp[, c("Model", "PercCov")]

# Percentages of genetic variance explained by the environmental covariables per management:
tmp <- aggregate(results, cbind(Covariables, LOF, Residual, Total) ~ Management, FUN = mean)
tmp$PercCov <- round(tmp$Covariables / (tmp$Covariables + tmp$LOF), 2)
tmp[, c("Management", "PercCov")]

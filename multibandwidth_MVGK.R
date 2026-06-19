library(asreml)
asreml.options(workspace = "5000mb", pworkspace = "5000mb", maxit = 50)
trace <- TRUE

# Seed:
set.seed(1994)
envs <- levels(readRDS("Briwecs/data/datalist.SE.rds")$ydata$Env)

# Loading kinship, env. correlation matrix and full data:
datalist <- readRDS("Briwecs/data/datalist.SE.rds")
d.full <- datalist$ydata

# Kinship and environmental matrices:
K <- readRDS("Briwecs/data/K.SE.rds")[levels(d.full$Gen), levels(d.full$Gen)]
EC <- datalist$EC[levels(d.full$Env), levels(d.full$Env)]
ED <- datalist$ED[levels(d.full$Env), levels(d.full$Env)]

vars.init.mv <- aggregate(d.full, GY ~ ManEnv, FUN = function(x) var(x) / 10)
vars.init.mv <- vars.init.mv$GY[match(levels(d.full$ManEnv), vars.init.mv$ManEnv)]
vars.init.sv <- rep(var(na.omit(d.full$GY)) / 20, 2)

## Multi-bandwidth, multi-var Gaussian kernel model ====
vf <- function(order, kappa) {
  # kappa[1] = variance of M1:E1
  # kappa[2] = variance of M1:E2
  #   ...
  # kappa[p * q] = variance of Mp:Eq
  # kappa[p * q + 1] = correlation between M1 and M2
  # kappa[p * q + 2] = bandwidth parameter of the Gaussian kernel for Env at M1
  # kappa[p * q + 3] = bandwidth parameter of the Gaussian kernel for Env at M2
  # kappa[p * q + 4] = bandwidth parameter of the Gaussian kernel for Env at M1:M2
  # Number of managements:
  n.mans <- order / nrow(ED)
  
  # The correlation matrix of the Man levels (specify manually!):
  Rm <- matrix(1, 2, 2)
  Rm[1, 2] <- Rm[2, 1] <- kappa[order + 1]
  
  # The full covariance matrix:
  S <- outer(sqrt(kappa[1:order]), sqrt(kappa[1:order]))
  # R <- kronecker(Rm, exp(-kappa[order + 2] * ED))
  R <- rbind(cbind(exp(-kappa[order + 2] * ED),               Rm[1, 2] * exp(-kappa[order + 4] * ED)),
             cbind(Rm[2, 1] * exp(-kappa[order + 4] * ED),    exp(-kappa[order + 3] * ED)))
  V <- S * R
  
  # Derivative wrt kappa[p * q + 1]
  # Indicator matrix of where kappa[p * q + 1] is present:
  I <- matrix(1, nrow(Rm), ncol(Rm))
  diag(I) = 0
  I <- kronecker(I, matrix(1, nrow(ED), ncol(ED)))
  dkrm <- (S * I) * kronecker(matrix(1, nrow(Rm), ncol(Rm)), exp(-kappa[order + 4] * ED))
  
  # Derivative wrt kappa[p * q + 2]
  dkh1 <- S * kronecker(matrix(c(1, 0, 0, 0), 2, 2), -ED * exp(-kappa[order + 2] * ED))
  
  # Derivative wrt kappa[p * q + 3]
  dkh2 <- S * kronecker(matrix(c(0, 0, 0, 1), 2, 2), -ED * exp(-kappa[order + 3] * ED))
  
  # Derivative wrt kappa[p * q + 4]
  dkh3 <- S * kronecker(matrix(c(0, kappa[order + 1], kappa[order + 1], 0), 2, 2), -ED * exp(-kappa[order + 4] * ED))
  
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
    deriv <- 0.5 * I * tmp * R
    deriv[dk, dk] <- 1
    varderivs[[dk]] <- deriv
  }
  # cat(kappa, "\n\n")
  return(c(list(V), varderivs, list(dkrm, dkh1, dkh2, dkh3)))
}

init <- c(vars.init.mv, 0.50, 0.10, 0.20, 0.15)
type <- c(rep("V", nrow(ED) * length(levels(d.full$Man))), "R", "V", "V", "V")
con <- c(rep("P", nrow(ED) * length(levels(d.full$Man))), "U", "P", "P", "P")
mod.mb.mvar.GK <- asreml(fixed = GY ~ -1 + ManEnv,
                         random = ~ own(ManEnv, "vf", init, type, con):vm(Gen, K),
                         residual = ~ units,
                         data = d.full,
                         trace = trace,
                         weights = weight,
                         family = asr_gaussian(dispersion = 1.0))

## Multi-bandwidth, single-var Gaussian kernel model ====
vf <- function(order, kappa) {
  # kappa[1] = variance of M1
  # kappa[2] = variance of M2
  # kappa[3] = correlation between M1 and M2
  # kappa[4] = bandwidth parameter of the Gaussian kernel for Env at M1
  # kappa[5] = bandwidth parameter of the Gaussian kernel for Env at M2
  # kappa[6] = bandwidth parameter of the Gaussian kernel for Env at M1:M2
  # The correlation matrix of the Man levels (specify manually!):
  Rm <- matrix(1, 2, 2)
  Rm[1, 2] <- Rm[2, 1] <- kappa[3]
  S <- outer(sqrt(kappa[1:ncol(Rm)]), sqrt(kappa[1:ncol(Rm)]))
  S <- kronecker(S, matrix(1, nrow(ED), ncol(ED)))
  
  # The full covariance matrix:
  # R <- kronecker(Rm, exp(-kappa[order + 2] * ED))
  R <- rbind(cbind(exp(-kappa[4] * ED),               Rm[1, 2] * exp(-kappa[6] * ED)),
             cbind(Rm[2, 1] * exp(-kappa[6] * ED),    exp(-kappa[5] * ED)))
  V <- S * R
  
  # Derivative wrt kappa[3]
  # Indicator matrix of where kappa[3] is present:
  I <- matrix(1, nrow(Rm), ncol(Rm))
  diag(I) = 0
  I <- kronecker(I, matrix(1, nrow(ED), ncol(ED)))
  dkrm <- (S * I) * kronecker(matrix(1, nrow(Rm), ncol(Rm)), exp(-kappa[6] * ED))
  
  # Derivative wrt kappa[4]
  dkh1 <- S * kronecker(matrix(c(1, 0, 0, 0), 2, 2), -ED * exp(-kappa[4] * ED))
  
  # Derivative wrt kappa[5]
  dkh2 <- S * kronecker(matrix(c(0, 0, 0, 1), 2, 2), -ED * exp(-kappa[5] * ED))
  
  # Derivative wrt kappa[6]
  dkh3 <- S * kronecker(matrix(c(0, kappa[3], kappa[3], 0), 2, 2), -ED * exp(-kappa[6] * ED))
  
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
  # cat(kappa, "\n\n")
  return(c(list(V), varderivs, list(dkrm, dkh1, dkh2, dkh3)))
}

init <- c(vars.init.sv, 0.50, 0.10, 0.20, 0.15)
type <- c("V", "V", "R", "V", "V", "V")
con <- c("P", "P", "U", "P", "P", "P")
mod.mb.svar.GK <- asreml(fixed = GY ~ -1 + ManEnv,
                         random = ~ own(ManEnv, "vf", init, type, con):vm(Gen, K),
                         residual = ~ units,
                         data = d.full,
                         trace = trace,
                         weights = weight,
                         family = asr_gaussian(dispersion = 1.0))


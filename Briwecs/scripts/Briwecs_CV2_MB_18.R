#!/usr/bin/env Rscript
# Prelims ====
# ASREML options:
library(asreml)
library(tictoc)
asreml.options(workspace = "5000mb", pworkspace = "5000mb", maxit = 25)
save.models <- FALSE
trace <- FALSE
arr.index <- 18
arr.dim <- 25

# Seed:
# set.seed(1994) # We don't need seeds because we are re-using the training and test sets that were randomly generated in the scripts for the other models.
# seeds <- floor(runif(arr.dim, 1000, 2000))
arr.slurm <- matrix(1:250, ncol = arr.dim)
envs <- levels(readRDS("Briwecs/data/datalist.SE.rds")$ydata$Env)

# Loading kinship, env. correlation matrix and full data:
datalist <- readRDS("Briwecs/data/datalist.SE.rds")
d.full <- datalist$ydata

# Kinship and environmental matrices:
K <- readRDS("Briwecs/data/K.SE.rds")[levels(d.full$Gen), levels(d.full$Gen)]
EC <- datalist$EC[levels(d.full$Env), levels(d.full$Env)]
ED <- datalist$ED[levels(d.full$Env), levels(d.full$Env)]

# Result storage:
runs <- nrow(arr.slurm) * ncol(arr.slurm) # Number of random training and test sets
reps <- 2 # Number of trials that each non-check genotype occurs in
models <- c("mb.svar.GK",             # multi-bandwidth single variance Gaussian kernel for E
            "mb.mvar.GK"              # multi-bandwidth multiple variance Gaussian kernel for E
)

n.checks <- (1:nrow(K))[((nrow(K) - (1:nrow(K)))/length(levels(d.full$Env))) %% 1 == 0][1:8]

results <- expand.grid(Run = 1:runs,
                       Checks = n.checks,
                       Model = models,
                       Man = levels(d.full$Man),
                       Env = rownames(ED),
                       cor_pearson = NA,     # Pearson correlation
                       cor_spearman = NA,    # Spearman rank correlation
                       RMSE = NA,            # Root mean squared error
                       MAE = NA,
                       Converged = NA)             # Median prediction standard error

comptimes <- expand.grid(Run = 1:runs,
                         Checks = n.checks,
                         Model = models,
                         Time = NA)

# Saving model fits in case we want to look at something:
if (save.models) {
  mod.fits <- vector("list", length(models))
  names(mod.fits) <- models
  for (i in 1:length(mod.fits)) {
    mod.fits[[i]] <- vector("list", runs)
  }
}

# Structure for saving the bandwidths:
bandwidths <- expand.grid(Run = 1:runs,
                          Checks = n.checks,
                          Model = models,
                          HN = NA,
                          LN = NA,
                          HN_LN = NA)

datasets <- arr.slurm[, arr.index]

cvsets <- readRDS(sprintf("Briwecs/results/CV2/cvsets.CV2.SE.%d-%d.rds", datasets[1], datasets[length(datasets)]))

nc <- n.checks[1]
run <- datasets[1]
start <- Sys.time()
for (nc in n.checks) {
  for (run in datasets) {
    # Setting the training and test sets that were already randomly generated for the previous models:
    train.set <- cvsets[[sprintf("nc%d", nc)]][[sprintf("run%d", run)]]$train.set
    test.set <- cvsets[[sprintf("nc%d", nc)]][[sprintf("run%d", run)]]$test.set
    
    # Making the datasets:
    d.train <- d.test <- d.full
    d.train[d.train$EnvGen %in% test.set, "GY"] <- NA
    d.test <- droplevels(d.full[d.full$EnvGen %in% test.set,])
    
    # Getting initial values for the MV and SV kernel model variances:
    # These match what asreml would use for corgh(ManEnv) and corgh(Man):corg(Env) models
    vars.init.mv <- aggregate(d.train, GY ~ ManEnv, FUN = function(x) var(x) / 10)
    vars.init.mv <- vars.init.mv$GY[match(levels(d.train$ManEnv), vars.init.mv$ManEnv)]
    vars.init.sv <- rep(var(na.omit(d.train$GY)) / 20, 2)
    
    # Models ====
    ## Multi-bandwidth, single-var Gaussian kernel model ====
    cat(sprintf("Fitting nc = %d run = %d, multi bandwidth, single var Gaussian kernel model...\n", nc, run))
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
        tmp <- kronecker(tmp, matrix(1, nrow(ED), ncol(ED)))
        deriv <- 0.5 * I * tmp * R
        # deriv <- 0.5 * I * kronecker(tmp * Rm, exp(-kappa[4] * ED))
        deriv[((dk - 1) * nrow(ED) + 1):(dk * nrow(ED)), ((dk - 1) * nrow(ED) + 1):(dk * nrow(ED))] <-
          deriv[((dk - 1) * nrow(ED) + 1):(dk * nrow(ED)), ((dk - 1) * nrow(ED) + 1):(dk * nrow(ED))] * 2
        varderivs[[dk]] <- deriv
      }
      # cat(kappa, "\n\n")
      return(c(list(V), varderivs, list(dkrm, dkh1, dkh2, dkh3)))
    }
    
    try({
      tic()
      init <- c(vars.init.sv, 0.50, 0.10, 0.10, 0.10)
      type <- c("V", "V", "R", "V", "V", "V")
      con <- c("P", "P", "U", "P", "P", "P")
      mod.svar.GK <- asreml(fixed = GY ~ -1 + ManEnv,
                            random = ~ own(ManEnv, "vf", init, type, con):vm(Gen, K),
                            residual = ~ units,
                            data = d.train,
                            trace = trace,
                            weights = weight,
                            family = asr_gaussian(dispersion = 1.0))
      
      pred.svar.GK <- as.data.frame(mod.svar.GK$coefficients$random)
      names(pred.svar.GK) <- "predicted.value"
      pred.svar.GK$Man <- as.factor(gsub("ManEnv_(.*):(.*):vm\\(Gen, K\\)_(.*)", "\\1", rownames(pred.svar.GK)))
      pred.svar.GK$Env <- as.factor(gsub("ManEnv_(.*):(.*):vm\\(Gen, K\\)_(.*)", "\\2", rownames(pred.svar.GK)))
      pred.svar.GK$Gen <- as.factor(gsub("ManEnv_(.*):(.*):vm\\(Gen, K\\)_(.*)", "\\3", rownames(pred.svar.GK)))
      pred.svar.GK$ManEnvGen <- as.factor(paste(pred.svar.GK$Man, pred.svar.GK$Env, pred.svar.GK$Gen, sep = ":"))
      pred.svar.GK <- pred.svar.GK[match(d.test$ManEnvGen, pred.svar.GK$ManEnvGen),]
      toc(log = TRUE)
      
      # Storing elapsed time:
      comptimes[comptimes$Run == run & comptimes$Model == "mb.svar.GK" & comptimes$Checks == nc,
                "Time"] <- as.numeric(unlist(lapply(tic.log(format = FALSE), function(x) x$toc - x$tic)))
      tic.clearlog()
      
      # Storing model fit:
      if (save.models) {
        mod.fits[["mb.svar.GK"]][[run]] <- mod.svar.GK
      }
      
      # Calculating and storing accuracies:
      for (man in levels(d.train$Man)) {
        for (e in levels(d.train$Env)) {
          preds <- droplevels(pred.svar.GK[pred.svar.GK$Env == e & pred.svar.GK$Man == man, c("Gen", "predicted.value")])
          target <- droplevels(d.test[d.test$Env == e & d.test$Man == man, c("Gen", "GY")])
          target <- target[match(preds$Gen, target$Gen),]
          
          # Accuracy measures:
          target$GY <- scale(target$GY, scale = FALSE)
          cor_pearson <- cor(preds$predicted.value, target$GY)
          cor_spearman <- cor(method = "spearman", preds$predicted.value, target$GY)
          RMSE <- sqrt(mean((preds$predicted.value - target$GY)^2))
          MAE <- mean(abs(preds$predicted.value - target$GY))
          
          result.row <- which(results$Run == run & results$Model == "mb.svar.GK" & results$Env == e & results$Man == man & results$Checks == nc)
          results[result.row, c("cor_pearson", "cor_spearman", "RMSE", "MAE", "Converged")] <-
            c(cor_pearson, cor_spearman, RMSE, MAE, mod.svar.GK$converge)
        }
      }
      
      # Storing bandwidths:
      bandwidths.row <- which(bandwidths$Run == run & bandwidths$Model == "mb.svar.GK" & bandwidths$Checks == nc)
      bandwidths[bandwidths.row, c("HN", "LN", "HN_LN")] <- summary(mod.svar.GK)$varcomp$component[4:6]
    })
    toc(); tic.clearlog()
    
    ## Multi-bandwidth, multi-var Gaussian kernel model ====
    cat(sprintf("Fitting nc = %d run = %d, multi bandwidth, multi var Gaussian kernel model...\n", nc, run))
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
    
    try({
      tic()
      init <- c(vars.init.mv, 0.50, 0.10, 0.10, 0.10)
      type <- c(rep("V", nrow(ED) * length(levels(d.full$Man))), "R", "V", "V", "V")
      con <- c(rep("P", nrow(ED) * length(levels(d.full$Man))), "U", "P", "P", "P")
      mod.mvar.GK <- asreml(fixed = GY ~ -1 + ManEnv,
                            random = ~ own(ManEnv, "vf", init, type, con):vm(Gen, K),
                            residual = ~ units,
                            data = d.train,
                            trace = trace,
                            weights = weight,
                            family = asr_gaussian(dispersion = 1.0))
      
      pred.mvar.GK <- as.data.frame(mod.mvar.GK$coefficients$random)
      names(pred.mvar.GK) <- "predicted.value"
      pred.mvar.GK$Man <- as.factor(gsub("ManEnv_(.*):(.*):vm\\(Gen, K\\)_(.*)", "\\1", rownames(pred.mvar.GK)))
      pred.mvar.GK$Env <- as.factor(gsub("ManEnv_(.*):(.*):vm\\(Gen, K\\)_(.*)", "\\2", rownames(pred.mvar.GK)))
      pred.mvar.GK$Gen <- as.factor(gsub("ManEnv_(.*):(.*):vm\\(Gen, K\\)_(.*)", "\\3", rownames(pred.mvar.GK)))
      pred.mvar.GK$ManEnvGen <- as.factor(paste(pred.mvar.GK$Man, pred.mvar.GK$Env, pred.mvar.GK$Gen, sep = ":"))
      pred.mvar.GK <- pred.mvar.GK[match(d.test$ManEnvGen, pred.mvar.GK$ManEnvGen),]
      toc(log = TRUE)
      
      # Storing elapsed time:
      comptimes[comptimes$Run == run & comptimes$Model == "mb.mvar.GK" & comptimes$Checks == nc,
                "Time"] <- as.numeric(unlist(lapply(tic.log(format = FALSE), function(x) x$toc - x$tic)))
      tic.clearlog()
      
      # Storing model fit:
      if (save.models) {
        mod.fits[["mb.mvar.GK"]][[run]] <- mod.mvar.GK
      }
      
      # Calculating and storing accuracies:
      for (man in levels(d.train$Man)) {
        for (e in levels(d.train$Env)) {
          preds <- droplevels(pred.mvar.GK[pred.mvar.GK$Env == e & pred.mvar.GK$Man == man, c("Gen", "predicted.value")])
          target <- droplevels(d.test[d.test$Env == e & d.test$Man == man, c("Gen", "GY")])
          target <- target[match(preds$Gen, target$Gen),]
          
          # Accuracy measures:
          target$GY <- scale(target$GY, scale = FALSE)
          cor_pearson <- cor(preds$predicted.value, target$GY)
          cor_spearman <- cor(method = "spearman", preds$predicted.value, target$GY)
          RMSE <- sqrt(mean((preds$predicted.value - target$GY)^2))
          MAE <- mean(abs(preds$predicted.value - target$GY))
          
          result.row <- which(results$Run == run & results$Model == "mb.mvar.GK" & results$Env == e & results$Man == man & results$Checks == nc)
          results[result.row, c("cor_pearson", "cor_spearman", "RMSE", "MAE", "Converged")] <-
            c(cor_pearson, cor_spearman, RMSE, MAE, mod.mvar.GK$converge)
        }
      }
      
      # Storing bandwidths:
      bandwidths.row <- which(bandwidths$Run == run & bandwidths$Model == "mb.mvar.GK" & bandwidths$Checks == nc)
      bandwidths[bandwidths.row, c("HN", "LN", "HN_LN")] <- summary(mod.mvar.GK)$varcomp$component[30:32]
    })
    toc(); tic.clearlog()
  }
}
end <- Sys.time()
end - start

saveRDS(results, sprintf("Briwecs/results/CV2/results.MB.CV2.SE.%d-%d.rds", datasets[1], datasets[length(datasets)]))
saveRDS(comptimes, sprintf("Briwecs/results/CV2/comptimes.MB.CV2.SE.%d-%d.rds", datasets[1], datasets[length(datasets)]))
saveRDS(bandwidths, sprintf("Briwecs/results/CV2/bandwidths.MB.CV2.SE.%d-%d.rds", datasets[1], datasets[length(datasets)]))
if (save.models) {
  saveRDS(mod.fits, sprintf("Briwecs/results/CV2/modfits.MB.CV2.SE.%d-%d.rds", datasets[1], datasets[length(datasets)]))
}






#!/usr/bin/env Rscript
# Prelims ====
# ASREML options:
library(asreml)
library(tictoc)
asreml.options(workspace = "5000mb", pworkspace = "5000mb", maxit = 25)
save.models <- FALSE
trace <- FALSE
arr.index <- 7
arr.dim <- 15

# Seed:
set.seed(1994)
seeds <- floor(runif(arr.dim, 1000, 2000))
arr.slurm <- matrix(1:150, ncol = arr.dim)
envs <- levels(readRDS("DROPS/data/ydata.rds")$Env)

# Loading kinship, env. correlation matrix and full data:
d.full <- readRDS("DROPS/data/ydata.rds")
d.full <- droplevels(d.full[d.full$Env %in% envs, c("Man", "Env", "Gen", "ManEnv", "ManGen", "EnvGen", "ManEnvGen", "GY")])

# Kinship and environmental matrices:
K <- readRDS("DROPS/data/K.rds")[levels(d.full$Gen), levels(d.full$Gen)]
EC <- readRDS("DROPS/data/EC.rds")[levels(d.full$Env), levels(d.full$Env)]
ED <- readRDS("DROPS/data/ED.rds")[levels(d.full$Env), levels(d.full$Env)]

# Result storage:
runs <- nrow(arr.slurm) * ncol(arr.slurm) # Number of random training and test sets
reps <- 2 # Number of trials that each non-check genotype occurs in
models <- c("ME",                  # No GxE, only a random main effect for genotype
            paste0("FA", c(1, 2, 3)), # Factor analytic models for the E term
            
            "svar.EC",             # single variance relmat for E
            "svar.GK",             # single variance Gaussian kernel for E
            
            "mvar.EC",             # multi variance relmat for E
            "mvar.GK"
)

# n.checks <- c(1, 8, 15, 22, 50, 78, 120)
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

datasets <- arr.slurm[, arr.index]
nc <- n.checks[1]
run <- datasets[1]
start <- Sys.time()

# Setting seed for each array job:
set.seed(seeds[arr.index])
for (nc in n.checks) {
  for (run in datasets) {
    # Generating data ====
    # Choose x random varieties as checks that occur in every trial:
    checks <- sample(levels(d.full$Gen), nc)
    gens <- setdiff(as.character(unique(d.full$Gen)), checks)
    train.set <- character(length(gens) * reps)
    
    while (train.set[length(train.set)] == "") {
      envs <- as.character(unique(d.full$Env)); names(envs) <- envs
      envs.count <- numeric(length(envs)); names(envs.count) <- envs
      i <- 1
      try({
        for (g in rep(setdiff(gens, checks), reps)) {
          e <- sample(envs, 1)
          q <- 1
          while (any(paste(e, g, sep = ":") %in% train.set) | any(envs.count[e] > (length(gens) * reps / nrow(ED) - 1))) {
            stopifnot("Sparse MET sampling got stuck... Trying again..." = q < 100)
            e <- sample(envs, 1)
            q <- q + 1
          }
          train.set[i] <- paste(e, g, sep = ":")
          i <- i + 1
          envs.count[e] <- envs.count[e] + 1
          envs.full <- names(envs.count)[which(envs.count == (length(gens) * reps / nrow(ED)))]
          envs <- envs[setdiff(envs, envs.full)]
        }
      })
    }
    checks <- paste(rep(as.character(unique(d.full$Env)), times = length(checks)), rep(checks, each = nrow(ED)), sep = ":")
    train.set <- c(train.set, checks)
    test.set <- setdiff(levels(d.full$EnvGen), train.set)
    
    # Making the datasets:
    d.train <- d.test <- d.full
    d.train[d.train$EnvGen %in% test.set, c("GY")] <- NA
    d.test <- droplevels(d.full[d.full$EnvGen %in% test.set,])
    
    # Models ====
    ##  ME model ====
    cat(sprintf("Fitting nc = %d run = %d, ME model...\n", nc, run))
    try({
      tic()
      mod.ME <- asreml(GY ~ -1 + ManEnv,
                       random = ~ vm(Gen, K),
                       residual = ~ units,
                       data = d.train,
                       trace = trace)
      
      pred.ME <- as.data.frame(mod.ME$coefficients$random)
      names(pred.ME) <- "predicted.value"
      pred.ME$Gen <- as.factor(gsub("vm\\(Gen, K\\)_(.*)", "\\1", rownames(pred.ME)))
      toc(log = TRUE)
      
      # Storing elapsed time:
      comptimes[comptimes$Run == run & comptimes$Model == "ME" & comptimes$Checks == nc,
                "Time"] <- as.numeric(unlist(lapply(tic.log(format = FALSE), function(x) x$toc - x$tic)))
      tic.clearlog()
      
      # Storing model fit:
      if (save.models) {
        mod.fits$ME[[run]] <- mod.ME
      }
      
      # Calculating and storing accuracies:
      for (man in levels(d.train$Man)) {
        for (e in levels(d.train$Env)) {
          
          preds <- pred.ME
          target <- droplevels(d.test[d.test$Env == e & d.test$Man == man, c("Gen", "GY")])
          preds <- preds[match(target$Gen, preds$Gen),]
          
          target$GY <- scale(target$GY, scale = FALSE)
          cor_pearson <- cor(preds$predicted.value, target$GY)
          cor_spearman <- cor(method = "spearman", preds$predicted.value, target$GY)
          RMSE <- sqrt(mean((preds$predicted.value - target$GY)^2))
          MAE <- mean(abs(preds$predicted.value - target$GY))
          
          result.row <- which(results$Run == run & results$Model == "ME" & results$Env == e & results$Man == man & results$Checks == nc)
          results[result.row, c("cor_pearson", "cor_spearman", "RMSE", "MAE", "Converged")] <-
            c(cor_pearson, cor_spearman, RMSE, MAE, mod.ME$converge)
        }
      }
    })
    # Extra toc in case the fit fails:
    toc(); tic.clearlog()
    
    
    
    ## FA models ====
    m <- 1
    for (m in c(1, 2, 3)) {
      cat(sprintf("Fitting nc = %d run = %d, FA%d model...\n", nc, run, m))
      try({
        tic()
        mod.FA <- asreml(GY ~ -1 + ManEnv,
                         random = ~ fa(ManEnv, m):vm(Gen, K),
                         residual = ~ units,
                         data = d.train,
                         trace = trace)
        
        pred.FA <- as.data.frame(mod.FA$coefficients$random[1:nrow(d.train),])
        names(pred.FA) <- "predicted.value"
        pred.FA$Man <- as.factor(gsub("fa\\(ManEnv, m\\)_(.*):(.*):vm\\(Gen, K\\)_(.*)", "\\1", rownames(pred.FA)))
        pred.FA$Env <- as.factor(gsub("fa\\(ManEnv, m\\)_(.*):(.*):vm\\(Gen, K\\)_(.*)", "\\2", rownames(pred.FA)))
        pred.FA$Gen <- as.factor(gsub("fa\\(ManEnv, m\\)_(.*):(.*):vm\\(Gen, K\\)_(.*)", "\\3", rownames(pred.FA)))
        pred.FA$ManEnvGen <- as.factor(paste(pred.FA$Man, pred.FA$Env, pred.FA$Gen, sep = ":"))
        pred.FA <- pred.FA[match(d.test$ManEnvGen, pred.FA$ManEnvGen),]
        toc(log = TRUE)
        
        # Storing elapsed time:
        comptimes[comptimes$Run == run & comptimes$Model == sprintf("FA%d", m) & comptimes$Checks == nc,
                  "Time"] <- as.numeric(unlist(lapply(tic.log(format = FALSE), function(x) x$toc - x$tic)))
        tic.clearlog()
        
        # Storing model fit:
        if (save.models) {
          mod.fits[[sprintf("FA%d", m)]][[run]] <- mod.FA
        }
        
        # Calculating and storing accuracies:
        for (man in levels(d.train$Man)) {
          for (e in levels(d.train$Env)) {
            
            preds <- droplevels(pred.FA[pred.FA$Env == e & pred.FA$Man == man, c("Gen", "predicted.value")])
            target <- droplevels(d.test[d.test$Env == e & pred.FA$Man == man, c("Gen", "GY")])
            target <- target[match(preds$Gen, target$Gen),]
            
            # Accuracy measures:
            target$GY <- scale(target$GY, scale = FALSE)
            cor_pearson <- cor(preds$predicted.value, target$GY)
            cor_spearman <- cor(method = "spearman", preds$predicted.value, target$GY)
            RMSE <- sqrt(mean((preds$predicted.value - target$GY)^2))
            MAE <- mean(abs(preds$predicted.value - target$GY))
            
            result.row <- which(results$Run == run & results$Model == sprintf("FA%d", m) & results$Env == e & results$Man == man & results$Checks == nc)
            results[result.row, c("cor_pearson", "cor_spearman", "RMSE", "MAE", "Converged")] <-
              c(cor_pearson, cor_spearman, RMSE, MAE, mod.FA$converge)
          }
        }
      })
      # Extra toc in case the fit fails:
      toc(); tic.clearlog()
    }
    
    
    
    ## Single-var relmat model ====
    cat(sprintf("Fitting nc = %d run = %d, single var relmat model...\n", nc, run))
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
      
      # cat(kappa, "\n\n")
      
      return(c(list(V),
               varderivs,
               list(dkcorr)))
    }
    
    try({
      tic()
      init <- c(0.1, 0.1, 0.1)
      type <- c("V", "V", "R")
      con <- c("P", "P", "U")
      mod.svar.EC <- asreml(fixed = GY ~ -1 + ManEnv,
                            random = ~ own(ManEnv, "vf", init, type, con):vm(Gen, K),
                            residual = ~ units,
                            data = d.train,
                            trace = trace)
      
      pred.svar.EC <- as.data.frame(mod.svar.EC$coefficients$random)
      names(pred.svar.EC) <- "predicted.value"
      pred.svar.EC$Man <- as.factor(gsub("ManEnv_(.*):(.*):vm\\(Gen, K\\)_(.*)", "\\1", rownames(pred.svar.EC)))
      pred.svar.EC$Env <- as.factor(gsub("ManEnv_(.*):(.*):vm\\(Gen, K\\)_(.*)", "\\2", rownames(pred.svar.EC)))
      pred.svar.EC$Gen <- as.factor(gsub("ManEnv_(.*):(.*):vm\\(Gen, K\\)_(.*)", "\\3", rownames(pred.svar.EC)))
      pred.svar.EC$ManEnvGen <- as.factor(paste(pred.svar.EC$Man, pred.svar.EC$Env, pred.svar.EC$Gen, sep = ":"))
      pred.svar.EC <- pred.svar.EC[match(d.test$ManEnvGen, pred.svar.EC$ManEnvGen),]
      toc(log = TRUE)
      
      # Storing elapsed time:
      comptimes[comptimes$Run == run & comptimes$Model == "svar.EC" & comptimes$Checks == nc,
                "Time"] <- as.numeric(unlist(lapply(tic.log(format = FALSE), function(x) x$toc - x$tic)))
      tic.clearlog()
      
      # Storing model fit:
      if (save.models) {
        mod.fits[["svar.EC"]][[run]] <- mod.svar.EC
      }
      
      # Calculating and storing accuracies:
      for (man in levels(d.train$Man)) {
        for (e in levels(d.train$Env)) {
          preds <- droplevels(pred.svar.EC[pred.svar.EC$Env == e & pred.svar.EC$Man == man, c("Gen", "predicted.value")])
          target <- droplevels(d.test[d.test$Env == e & d.test$Man == man, c("Gen", "GY")])
          target <- target[match(preds$Gen, target$Gen),]
          
          # Accuracy measures:
          target$GY <- scale(target$GY, scale = FALSE)
          cor_pearson <- cor(preds$predicted.value, target$GY)
          cor_spearman <- cor(method = "spearman", preds$predicted.value, target$GY)
          RMSE <- sqrt(mean((preds$predicted.value - target$GY)^2))
          MAE <- mean(abs(preds$predicted.value - target$GY))
          
          result.row <- which(results$Run == run & results$Model == "svar.EC" & results$Env == e & results$Man == man & results$Checks == nc)
          results[result.row, c("cor_pearson", "cor_spearman", "RMSE", "MAE", "Converged")] <-
            c(cor_pearson, cor_spearman, RMSE, MAE, mod.svar.EC$converge)
        }
      }
    })
    toc(); tic.clearlog()
    
    
    
    ## Single-var Gaussian kernel model ====
    cat(sprintf("Fitting nc = %d run = %d, single var Gaussian kernel model...\n", nc, run))
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
      
      return(c(list(V),
               varderivs,
               list(dkcorr, dkbw)))
    }
    
    try({
      tic()
      init <- c(0.1, 0.1, 0.1, 0.1)
      type <- c("V", "V", "R", "V")
      con <- c("P", "P", "U", "P")
      mod.svar.GK <- asreml(fixed = GY ~ -1 + ManEnv,
                            random = ~ own(ManEnv, "vf", init, type, con):vm(Gen, K),
                            residual = ~ units,
                            data = d.train,
                            trace = trace)
      
      pred.svar.GK <- as.data.frame(mod.svar.GK$coefficients$random)
      names(pred.svar.GK) <- "predicted.value"
      pred.svar.GK$Man <- as.factor(gsub("ManEnv_(.*):(.*):vm\\(Gen, K\\)_(.*)", "\\1", rownames(pred.svar.GK)))
      pred.svar.GK$Env <- as.factor(gsub("ManEnv_(.*):(.*):vm\\(Gen, K\\)_(.*)", "\\2", rownames(pred.svar.GK)))
      pred.svar.GK$Gen <- as.factor(gsub("ManEnv_(.*):(.*):vm\\(Gen, K\\)_(.*)", "\\3", rownames(pred.svar.GK)))
      pred.svar.GK$ManEnvGen <- as.factor(paste(pred.svar.GK$Man, pred.svar.GK$Env, pred.svar.GK$Gen, sep = ":"))
      pred.svar.GK <- pred.svar.GK[match(d.test$ManEnvGen, pred.svar.GK$ManEnvGen),]
      toc(log = TRUE)
      
      # Storing elapsed time:
      comptimes[comptimes$Run == run & comptimes$Model == "svar.GK" & comptimes$Checks == nc,
                "Time"] <- as.numeric(unlist(lapply(tic.log(format = FALSE), function(x) x$toc - x$tic)))
      tic.clearlog()
      
      # Storing model fit:
      if (save.models) {
        mod.fits[["svar.GK"]][[run]] <- mod.svar.GK
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
          
          result.row <- which(results$Run == run & results$Model == "svar.GK" & results$Env == e & results$Man == man & results$Checks == nc)
          results[result.row, c("cor_pearson", "cor_spearman", "RMSE", "MAE", "Converged")] <-
            c(cor_pearson, cor_spearman, RMSE, MAE, mod.svar.GK$converge)
        }
      }
    })
    toc(); tic.clearlog()
    
    
    
    ## Multi-var relmat model ====
    cat(sprintf("Fitting nc = %d run = %d, multi var relmat model...\n", nc, run))
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
    
    try({
      tic()
      init <- c(rep(0.1, nrow(EC) * length(levels(d.train$Man))), 0.1)
      type <- c(rep("V", nrow(EC) * length(levels(d.train$Man))), "R")
      con <- c(rep("P", nrow(EC) * length(levels(d.train$Man))), "U")
      mod.mvar.EC <- asreml(fixed = GY ~ -1 + ManEnv,
                            random = ~ own(ManEnv, "vf", init, type, con):vm(Gen, K),
                            residual = ~ units,
                            data = d.train,
                            trace = trace)
      
      pred.mvar.EC <- as.data.frame(mod.mvar.EC$coefficients$random)
      names(pred.mvar.EC) <- "predicted.value"
      pred.mvar.EC$Man <- as.factor(gsub("ManEnv_(.*):(.*):vm\\(Gen, K\\)_(.*)", "\\1", rownames(pred.mvar.EC)))
      pred.mvar.EC$Env <- as.factor(gsub("ManEnv_(.*):(.*):vm\\(Gen, K\\)_(.*)", "\\2", rownames(pred.mvar.EC)))
      pred.mvar.EC$Gen <- as.factor(gsub("ManEnv_(.*):(.*):vm\\(Gen, K\\)_(.*)", "\\3", rownames(pred.mvar.EC)))
      pred.mvar.EC$ManEnvGen <- as.factor(paste(pred.mvar.EC$Man, pred.mvar.EC$Env, pred.mvar.EC$Gen, sep = ":"))
      pred.mvar.EC <- pred.mvar.EC[match(d.test$ManEnvGen, pred.mvar.EC$ManEnvGen),]
      toc(log = TRUE)
      
      # Storing elapsed time:
      comptimes[comptimes$Run == run & comptimes$Model == "mvar.EC" & comptimes$Checks == nc,
                "Time"] <- as.numeric(unlist(lapply(tic.log(format = FALSE), function(x) x$toc - x$tic)))
      tic.clearlog()
      
      # Storing model fit:
      if (save.models) {
        mod.fits[["mvar.EC"]][[run]] <- mod.mvar.EC
      }
      
      # Calculating and storing accuracies:
      for (man in levels(d.train$Man)) {
        for (e in levels(d.train$Env)) {
          preds <- droplevels(pred.mvar.EC[pred.mvar.EC$Env == e & pred.mvar.EC$Man == man, c("Gen", "predicted.value")])
          target <- droplevels(d.test[d.test$Env == e & d.test$Man == man, c("Gen", "GY")])
          target <- target[match(preds$Gen, target$Gen),]
          
          # Accuracy measures:
          target$GY <- scale(target$GY, scale = FALSE)
          cor_pearson <- cor(preds$predicted.value, target$GY)
          cor_spearman <- cor(method = "spearman", preds$predicted.value, target$GY)
          RMSE <- sqrt(mean((preds$predicted.value - target$GY)^2))
          MAE <- mean(abs(preds$predicted.value - target$GY))
          
          result.row <- which(results$Run == run & results$Model == "mvar.EC" & results$Env == e & results$Man == man & results$Checks == nc)
          results[result.row, c("cor_pearson", "cor_spearman", "RMSE", "MAE", "Converged")] <-
            c(cor_pearson, cor_spearman, RMSE, MAE, mod.mvar.EC$converge)
        }
      }
    })
    toc(); tic.clearlog()
    
    
    
    ## Multi-var Gaussian kernel model ====
    cat(sprintf("Fitting nc = %d run = %d, multi var Gaussian kernel model...\n", nc, run))
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
    
    try({
      tic()
      init <- c(rep(0.1, nrow(ED) * length(levels(d.train$Man))), 0.1, 0.1)
      type <- c(rep("V", nrow(ED) * length(levels(d.train$Man))), "R", "V")
      con <- c(rep("P", nrow(ED) * length(levels(d.train$Man))), "U", "P")
      mod.mvar.GK <- asreml(fixed = GY ~ -1 + ManEnv,
                            random = ~ own(ManEnv, "vf", init, type, con):vm(Gen, K),
                            residual = ~ units,
                            data = d.train,
                            trace = trace)
      
      pred.mvar.GK <- as.data.frame(mod.mvar.GK$coefficients$random)
      names(pred.mvar.GK) <- "predicted.value"
      pred.mvar.GK$Man <- as.factor(gsub("ManEnv_(.*):(.*):vm\\(Gen, K\\)_(.*)", "\\1", rownames(pred.mvar.GK)))
      pred.mvar.GK$Env <- as.factor(gsub("ManEnv_(.*):(.*):vm\\(Gen, K\\)_(.*)", "\\2", rownames(pred.mvar.GK)))
      pred.mvar.GK$Gen <- as.factor(gsub("ManEnv_(.*):(.*):vm\\(Gen, K\\)_(.*)", "\\3", rownames(pred.mvar.GK)))
      pred.mvar.GK$ManEnvGen <- as.factor(paste(pred.mvar.GK$Man, pred.mvar.GK$Env, pred.mvar.GK$Gen, sep = ":"))
      pred.mvar.GK <- pred.mvar.GK[match(d.test$ManEnvGen, pred.mvar.GK$ManEnvGen),]
      toc(log = TRUE)
      
      # Storing elapsed time:
      comptimes[comptimes$Run == run & comptimes$Model == "mvar.GK" & comptimes$Checks == nc,
                "Time"] <- as.numeric(unlist(lapply(tic.log(format = FALSE), function(x) x$toc - x$tic)))
      tic.clearlog()
      
      # Storing model fit:
      if (save.models) {
        mod.fits[["mvar.GK"]][[run]] <- mod.mvar.GK
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
          
          result.row <- which(results$Run == run & results$Model == "mvar.GK" & results$Env == e & results$Man == man & results$Checks == nc)
          results[result.row, c("cor_pearson", "cor_spearman", "RMSE", "MAE", "Converged")] <-
            c(cor_pearson, cor_spearman, RMSE, MAE, mod.mvar.GK$converge)
        }
      }
    })
    toc(); tic.clearlog()
  }
}
end <- Sys.time()
end - start

saveRDS(results, sprintf("DROPS/results/CV2/results.CV2.%d-%d.rds", datasets[1], datasets[length(datasets)]))
saveRDS(comptimes, sprintf("DROPS/results/CV2/comptimes.CV2.%d-%d.rds", datasets[1], datasets[length(datasets)]))
if (save.models) {
  saveRDS(mod.fits, sprintf("DROPS/results/CV2/modfits.CV2.%d-%d.rds", datasets[1], datasets[length(datasets)]))
}



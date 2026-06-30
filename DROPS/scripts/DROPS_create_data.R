# We will use environmental covariables to obtain a correlation and distance
# matrix of 15 rain-fed trials through Europe.
set.seed(1997)

# Creating or loading kinship:
if (!("K.rds" %in% list.files("DROPS/data/"))) {
  geno <- read.csv("DROPS/raw_data/7a-Genotyping_50K_41722.csv")
  rownames(geno) <- geno$Ind
  geno <- geno[, -1]
  gdata <- statgenGWAS::createGData(geno = geno)
  gdata <- statgenGWAS::codeMarkers(gdata, removeDuplicates = T, MAF = 0.05)
  K <- statgenGWAS::kinship(gdata$markers, "IBS")
  saveRDS(K, "DROPS/data/K.rds")
} else {
  K <- readRDS("DROPS/data/K.rds")
}

# Loading plot-level data:
ydata <- read.csv("DROPS/raw_data/2a-GrainYield_components_Plot_level.csv")

# Subsetting to relevant columns:
ydata <- ydata[, c("Experiment", "treatment", "Replicate", "block", "Row", "Column", "Variety_ID", "grain.yield", "type")]
ydata$Man <- ifelse(ydata$treatment == "rainfed", "R", "W")
ydata$Env <- substr(ydata$Experiment, 1, 5)
ydata$Gen <- ydata$Variety_ID
ydata <- droplevels(na.omit(ydata))
ydata <- droplevels(ydata[ydata$type == "CrossValidation",])
# Discard Deb13 because it only has the R management:
ydata <- droplevels(ydata[ydata$Env != "Deb13",])
ydata <- ydata[, c("Man", "Env", "Gen", "Replicate", "block", "Row", "Column", "grain.yield")]
colnames(ydata) <- c("Man", "Env", "Gen", "Rep", "Block", "Row", "Col", "GY")
ydata$ManEnv <- paste(ydata$Man, ydata$Env, sep = ":")

# Fitting single-trial models:
e <- unique(ydata$ManEnv)[1]
i <- 1
# ggplot(droplevels(ydata[ydata$ManEnv == e & ydata$Block %in% 1:8,]), aes(x = Col, y = Row, color = as.factor(Block), shape = as.factor(Rep))) +
#   geom_point(size = 3) + theme_classic()
if (!(all(c("d.corr.SE.rds", "K.corr.SE.rds") %in% list.files("DROPS/data")))) {
  for (e in unique(ydata$ManEnv)) {
    ydata.ss <- droplevels(ydata[ydata$ManEnv == e,])
    ydata.ss$R <- as.factor(ydata.ss$Row)
    ydata.ss$C <- as.factor(ydata.ss$Col)
    ydata.ss$Gen <- as.factor(ydata.ss$Gen)
    ydata.ss$Rep <- as.factor(ydata.ss$Rep)
    fit <- LMMsolve(fixed = GY ~ -1 + Gen,
                    random = ~ Rep + R + C,
                    spline = ~ spl2D(Row, Col, nseg = c(25, 70)),
                    data = ydata.ss, maxit = 250, trace = TRUE, tolerance = 1e-6)
    
    # Getting standard errors
    C <- as.matrix(fit$C) # Coefficient matrix
    X <- as.matrix(fit$X) # Fixed effects design matrix
    Z <- as.matrix(fit$Z) # Random effects design matrix
    nX <- ncol(as.matrix(fit$X)) # Number of fixed effects (BLUEs for management-genotype combinations)
    Czz <- C[(nX + 1):nrow(C), (nX + 1):ncol(C)] # Block of C corresponding to the random effects (Zt %*% Ri %*% Z + Gi)
    Ri <- as.matrix(fit$lRinv$residual) * (1 / fit$VarDf$Variance[fit$VarDf$VarComp == "residual"]) # Inverse residual covariance matrix
    Gi <- Czz - t(Z) %*% Ri %*% Z # Inverse random effects covariance matrix, see comment on Czz line
    V <- Z %*% solve(Gi) %*% t(Z) + solve(Ri) # Total covariance matrix
    Vi <- solve(V)
    omega <- solve(t(X) %*% Vi %*% X)
    omegai <- solve(omega)
    weights <- diag(omegai)[1:(nrow(omegai) - 3)]
    
    estimates <- coef(fit)
    BLUEs <- data.frame(Gen = gsub("Gen_(.*)", "\\1", names(estimates$Gen)),
                        GY = as.numeric(estimates$Gen))
    BLUEs$ManEnv <- e
    BLUEs$Man <- gsub("(.):.*", "\\1", e)
    BLUEs$Env <- gsub(".:(.*)", "\\1", e)
    BLUEs <- BLUEs[-1, ] # Discarding the first row because we have no intercept.
    if (length(weights) != nrow(BLUEs)) {
      stop("Something went wrong!")
    }
    BLUEs$weight <- weights
    reps <- as.data.frame(table(ydata.ss$Gen))
    BLUEs$reps <- reps[match(BLUEs$Gen, as.character(reps$Var1)), "Freq"]
    if (i == 1) {
      d.corr <- BLUEs[, c("Man", "Env", "Gen", "GY", "weight", "reps")]
    } else {
      d.corr <- rbind(d.corr, BLUEs[, c("Man", "Env", "Gen", "GY", "weight", "reps")])
    }
    i <- i + 1
  }
  d.corr <- d.corr[which(d.corr$Gen %in% rownames(K)),]
  d.corr$Man <- as.factor(d.corr$Man)
  d.corr$Env <- as.factor(d.corr$Env)
  d.corr$Gen <- as.factor(d.corr$Gen)
  d.corr$ManEnv <- as.factor(paste(as.character(d.corr$Man), as.character(d.corr$Env), sep = ":"))
  saveRDS(d.corr, "DROPS/data/d.corr.SE.rds")
  saveRDS(K, "DROPS/data/K.corr.SE.rds")
} else {
  d.corr <- readRDS("DROPS/data/d.corr.SE.rds")
  K <- readRDS("DROPS/data/K.corr.SE.rds")
}

# Now use d.corr instead of loading in the published BLUEs and continue as before:
ydata <- d.corr

# Loading and subsetting yield data:
# ydata <- read.csv("DROPS/raw_data/2b-GrainYield_components_BLUEs_level-1.csv")
# ydata <- ydata[which(ydata$Variety_ID %in% rownames(K)),]

# There is one trial that only has 238 genotypes (W:Kar13):
# aggregate(ydata, GY ~ ManEnv, FUN = length)

# We need to figure out which genotypes are present in all 28 trials:
present <- vector("list", length(levels(ydata$ManEnv)))
for (i in 1:length(present)) {
  names(present)[i] <- levels(ydata$ManEnv)[i]
  present[[i]] <- levels(droplevels(ydata[ydata$ManEnv == names(present)[i],])$Gen)
}
keep <- Reduce(intersect, present) # 234 genotypes that are present everywhere, let's just discard the 12 other genotypes.
ydata <- droplevels(ydata[ydata$Gen %in% keep,])

# Loading environmental data:
edata <- read.csv("DROPS/raw_data/1-Env_variables_daily.csv")

# Which loc-years do we have in ydata (14 loc-years)?
y.locyears <- unique(as.character(ydata$Env))

# We have edata for all 14 loc-years in ydata!
length(y.locyears[which(y.locyears %in% unique(edata$Env))])

# Subsetting edata:
edata <- edata[which(edata$Env %in% y.locyears),]

# Converting the dates to day of year:
edata$DOY <- as.numeric(strftime(edata$Date, format = "%j"))
edata$day <- numeric(nrow(edata))

# Checking whether we have data for all days from start to end:
for (i in setdiff(y.locyears, "Gra13")) {
  edata[which(edata$Env == i), "day"] <-
    edata[which(edata$Env == i), "DOY"] - min(edata[which(edata$Env == i), "DOY"]) + 1
  
  if (all.equal(edata[which(edata$Env == i), "day"], 1:length(edata[which(edata$Env == i), "DOY"]))) {
    cat(sprintf("%s is correct!\n", i))
  } else {
    cat(sprintf("%s is missing a date!\n", i))
  }
}

# Gra13 is annoying because it's Chile where sowing is done in November:
# So we check manually:
all.equal(c(327:(327+38), 1:104), edata[which(edata$Env == "Gra13"), "DOY"])
edata[which(edata$Env == "Gra13"), "day"] <- 1:length(edata[which(edata$Env == "Gra13"), "day"])

# Now we calculate corn growing degree days:
edata$GDD <- edata$cumulGDD <- numeric(nrow(edata))
gdd <- function(tmin, tmax) {
  tmin <- tmin * 1.8 + 32 # Convert to F
  tmax <- tmax * 1.8 + 32
  if (tmin < 50) tmin <- 50 # Minimum and maximum growing temperatures for corn
  if (tmax < 50) tmax <- 50
  if (tmin > 86) tmin <- 86
  if (tmax > 86) tmax <- 86
  GDD <- ((tmin + tmax) / 2) - 50
}

# Calculate GDDs for individual dates:
for (i in 1:nrow(edata)) {
  edata$GDD[i] <- gdd(tmin = edata$Tmin.air[i], tmax = edata$Tmax.air[i])
}

# Calculate cumulative GDDs:
for (env in unique(edata$Env)) {
  for (d in 1:nrow(edata[which(edata$Env == env),])) {
    edata[which(edata$Env == env & edata$day == d), "cumulGDD"] <-
      sum(edata[which(edata$Env == env & edata$day %in% 1:d), "GDD"])
  }
}

# Printing end of season GDDs (quite a bit of variation which is not ideal):
for (env in unique(edata$Env)) {
  cat(sprintf("Environment: %s\tEOF GDD: %s\n", env, round(max(edata[which(edata$Env == env), "cumulGDD"]))))
}

# Variables to potentially use for the environmental correlation matrix:
vars <- c("RHmin.air", "RHmax.air", "RHmean.air", "Raincum", "Windspeedmax",
          "ET0.air", "Rad", "Ri", "VPD.air", "VPD.apex", "Tmax.apex", "Tnight")

edata <- edata[, c("Env", "Date", "DOY", "day", "GDD", "cumulGDD", vars)]

# Cumulative GDD windows:
wins <- c("(0,100]", "(100,200]", "(200,300]", "(300,400]", "(400,500]", "(500,600]", "(600,700]", "(700,800]",
          "(800,900]", "(900,1e+03]", "(1e+03,1.1e+03]", "(1.1e+03,1.2e+03]", "(1.2e+03,1.3e+03]",
          "(1.3e+03,1.4e+03]", "(1.4e+03,1.5e+03]", "(1.5e+03,1.6e+03]", "(1.6e+03,1.7e+03]", "(1.7e+03,1.8e+03]",
          "(1.8e+03,1.9e+03]", "(1.9e+03,2e+03]", "(2e+03,2.1e+03]", "(2.1e+03,2.2e+03]", "(2.2e+03,2.3e+03]",
          "(2.3e+03,2.4e+03]", "(2.4e+03,2.5e+03]", "(2.5e+03,2.6e+03]", "(2.6e+03,2.7e+03]", "(2.7e+03,2.8e+03]",
          "(2.8e+03,2.9e+03]", "(2.9e+03,3e+03]", "(3e+03,3.1e+03]")

edata.fit <- expand.grid(env = unique(edata$Env),
                         var = vars,
                         window = wins,
                         AUC = NA,
                         stringsAsFactors = F)

# Scaling environmental variables:
edata[, 7:ncol(edata)] <- scale(edata[, 7:ncol(edata)])

# Fitting piecewiese constant regressions for each variable with breakpoints at each
# transition from one window to the next.
for (e in unique(edata.fit$env)) {
  for (t in unique(edata.fit$var)) {
    
    edata.ss <- na.omit(edata[which(edata$Env == e), c(t, "cumulGDD")])
    edata.ss$window <- cut(edata.ss$cumulGDD, breaks = seq(0, 3100, 100))
    fml <- paste0(t, " ~ window")
    fit <- lm(as.formula(fml), data = edata.ss)
    windows <- unique(as.character(edata.ss$window))
    intercepts <- coef(fit)
    intercepts[2:length(intercepts)] <- intercepts[2:length(intercepts)] + intercepts[1]
    names(intercepts) <- windows
    for (s in windows) {
      edata.fit[which(edata.fit$env == e & edata.fit$var == t & edata.fit$window == s), "AUC"] <-
        intercepts[s]
    }
  }
}

# Using only AUCs for stages that are present in all trials:
discard <- unique(edata.fit[which(is.na(edata.fit$AUC)), "window"])
edata.fit <- droplevels(edata.fit[which(!(edata.fit$window %in% discard)),])
any(is.na(edata.fit)) # FALSE
edata.fit.wide <- as.data.frame(tidyr::pivot_wider(edata.fit, names_from = 2:3, values_from = 4))
rownames(edata.fit.wide) <- edata.fit.wide$env

# Making and saving environmental correlation matrix:
X <- scale(t(as.matrix(edata.fit.wide[, -1])))
C <- (t(X) %*% X) / (ncol(edata.fit.wide[, -1]) - 1)
saveRDS(C, "DROPS/data/EC.SE.rds")

# Distance matrix:
X <- scale(as.matrix(edata.fit.wide[, -1]))
ED <- (as.matrix(dist(X, method = "euclidian"))^2) / ncol(X)
saveRDS(ED, "DROPS/data/ED.SE.rds")

# Making 100% sure the ordering is right:
K <- K[levels(ydata$Gen), levels(ydata$Gen)]
M <- matrix(0, 2, 2)
rownames(M) <- colnames(M) <- c("R", "W")
order <- rownames(kronecker(kronecker(M, ED, make.dimnames = TRUE), K, make.dimnames = TRUE))
ydata$ManEnvGen <- paste(ydata$Man, ydata$Env, ydata$Gen, sep = ":")
ydata$ManEnv <- paste(ydata$Man, ydata$Env, sep = ":")
ydata$ManGen <- paste(ydata$Man, ydata$Gen, sep = ":")
ydata$EnvGen <- paste(ydata$Env, ydata$Gen, sep = ":")
ydata <- ydata[match(order, ydata$ManEnvGen),]

ydata$Env <- factor(ydata$Env, levels = ydata$Env, labels = ydata$Env)
ydata$Man <- factor(ydata$Man, levels = ydata$Man, labels = ydata$Man)
ydata$Gen <- factor(ydata$Gen, levels = ydata$Gen, labels = ydata$Gen)
ydata$ManEnvGen <- factor(ydata$ManEnvGen, levels = ydata$ManEnvGen, labels = ydata$ManEnvGen)
ydata$ManEnv <- factor(ydata$ManEnv, levels = ydata$ManEnv, labels = ydata$ManEnv)
ydata$ManGen <- factor(ydata$ManGen, levels = ydata$ManGen, labels = ydata$ManGen)
ydata$EnvGen <- factor(ydata$EnvGen, levels = ydata$EnvGen, labels = ydata$EnvGen)

all(levels(ydata$Env) == rownames(ED))
all(levels(ydata$Gen) == rownames(K))
all(levels(ydata$ManEnv) == rownames(kronecker(M, ED, make.dimnames = TRUE)))
all(levels(ydata$ManEnvGen) == rownames(kronecker(kronecker(M, ED, make.dimnames = TRUE), K, make.dimnames = TRUE)))

all(levels(ydata$Man) == unique(as.character(ydata$Man)))
all(levels(ydata$Env) == unique(as.character(ydata$Env)))
all(levels(ydata$Gen) == unique(as.character(ydata$Gen)))

all(levels(ydata$ManEnvGen) == unique(as.character(ydata$ManEnvGen)))

all(levels(ydata$ManEnv) == unique(as.character(ydata$ManEnv)))
all(levels(ydata$ManGen) == unique(as.character(ydata$ManGen)))
all(levels(ydata$EnvGen) == unique(as.character(ydata$EnvGen)))

# Saving:
saveRDS(ydata, "DROPS/data/ydata.SE.rds")
saveRDS(K, "DROPS/data/K.SE.rds")



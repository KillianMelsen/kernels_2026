library(LMMsolver)
library(ggplot2)
library(asreml)
asreml.options(workspace = "5000mb", pworkspace = "5000mb", maxit = 50)
set.seed(1997)

# Cultivar info:
tmp <- read.csv("Briwecs/raw_data/BRIWECs_cultivar_info.csv")[, c("BRISONr", "genotype")]

# Some genotypes have annoying names so we maually change those:
tmp$genotype[which(tmp$genotype == "Triple Dirk \"S\"")] <- "Triple Dirk S"
tmp$genotype[which(tmp$genotype == "G\xf6tz")] <- "Gotz"
tmp$genotype[which(tmp$genotype == "T\xfcrkis")] <- "Turkis"

# Phenotypic data:
d <- read.csv("Briwecs/raw_data/BRIWECS_data_publication.csv", sep = ";")
d <- droplevels(d[d$BRISONr != "BRISONr_229",]) # Dropping this genotype
for (i in 1:nrow(d)) {
  if (grepl("BRISONr_.*", d$BRISONr[i])) {
    d$BRISONr[i] <- tmp$genotype[tmp$BRISONr == d$BRISONr[i]]
  }
}
rm(tmp, i)
d$BRISONr <- gsub(" ", "_", d$BRISONr)

# For pre-processing, we need phase 1 (2015-2017) and the RF managements (low/high nitrogen, with/without fungicide, rainfed):
# d <- droplevels(d[d$Treatment %in% c("LN_NF_RF", "HN_NF_RF") & d$Year %in% c(2015, 2016, 2017),])
d <- droplevels(d[d$Year %in% c(2015, 2016, 2017) & d$Treatment %in% c("LN_NF_RF", "HN_NF_RF", "LN_WF_RF", "HN_WF_RF"),])
d$Env <- as.factor(paste(d$Location, d$Year, sep = "_"))
d.full <- d
d <- d[, c("BRISONr", "Treatment", "Block", "Row", "Column", "Sowing_date", "Emergence_date", "BBCH87", "Seedyield", "Env")]
colnames(d) <- c("Gen", "Man", "Block", "Row", "Col", "Sowing", "Emergence", "Maturity", "GY", "Env")
d <- d[, c("Man", "Env", "Gen", "Block", "Row", "Col", "Sowing", "Emergence", "Maturity", "GY")]
# d$Man <- ifelse(d$Man == "LN_NF_RF", "LN", "HN")
d$Man <- as.factor(d$Man)
d$Gen <- as.factor(d$Gen)
d$ManEnv <- as.factor(paste(d$Man, d$Env, sep = ":"))
d$Block <- as.factor(d$Block)
d$RowCol <- as.factor(paste(d$Row, d$Col, sep = ":"))

# We will discard the RHH_2016 environment because there are multiple copies of
# row-col combinations in this trial...
# d <- droplevels(d[d$Env != "RHH_2016",])

# We will also discard RHH_2017 because it seems to only have 50 genotypes:
d <- droplevels(d[d$Env != "RHH_2017",])

# We will also discard GGE_2015 because it only has a single management:
d <- droplevels(d[d$Env != "GGE_2015",])

# We will also discard GGE_2016 and GGE_2017 because they seem to have 4 replicates, but only 2 blocks.
# Not entirely clear what happened there, because the paper only mentions that there should be two
# replicates. They "fixed" this, but now the GGE environments don't have the managements
# mentioned in the paper anymore...
d <- droplevels(d[d$Env != "GGE_2016",])
d <- droplevels(d[d$Env != "GGE_2017",])

# Some plotting to understand the trial layout:
# levels(d$Env)
# ggplot(droplevels(d[d$Env == "KIE_2017",]), aes(x = Col, y = Row, color = Man, shape = Block)) +
#   geom_point(size = 3) + theme_classic()
if (!("d.corr.rds" %in% list.files("Briwecs/data/"))) {
  BLUEs <- vector("list", length(levels(d$Env)))
  i <- 14
  for (i in 1:length(levels(d$Env))) {
    e <- levels(d$Env)[i]
    dss <- droplevels(d[d$Env == e,])
    dss$R <- as.factor(dss$Row)
    dss$C <- as.factor(dss$Col)
    dss$ManGen <- as.factor(paste(dss$Man, dss$Gen, sep = ":"))
    dss <- droplevels(dss[!(is.na(dss$GY)),])
    fit <- LMMsolve(fixed = GY ~ ManGen,
                    random = ~ Man:Block + R + C,
                    spline = ~ spl2D(Row, Col, nseg = c(25, 70)),
                    data = dss, maxit = 250, trace = TRUE, tolerance = 1e-6)
    
    pred <- obtainSmoothTrend(fit, newdata = dss, includeIntercept = TRUE)
  
    ggplot(pred, aes(x = Col, y = Row, fill = ypred)) +
      geom_tile(show.legend = TRUE) +
      scale_fill_gradientn(colours = topo.colors(100)) +
      coord_fixed() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    
    estimates <- coef(fit)
    BLUEs <- data.frame(ManGen = names(estimates$ManGen),
                        GY = as.numeric(estimates$ManGen))
    BLUEs$GY <- BLUEs$GY + as.numeric(estimates$`(Intercept)`)
    BLUEs$Man <- gsub("ManGen_(.*):(.*)", "\\1", BLUEs$ManGen)
    BLUEs$Gen <- gsub("ManGen_(.*):(.*)", "\\2", BLUEs$ManGen)
    BLUEs$Env <- e
    if (i == 1) {
      d.corr <- BLUEs[, c("Man", "Env", "Gen", "GY")]
    } else {
      d.corr <- rbind(d.corr, BLUEs[, c("Man", "Env", "Gen", "GY")])
    }
  }
  
  d.corr$Man <- as.factor(d.corr$Man)
  d.corr$Env <- as.factor(d.corr$Env)
  d.corr$Gen <- as.factor(d.corr$Gen)
  
  # Getting the weights (as the number of replicates):
  d.reps <- droplevels(d[!(is.na(d$GY)),])
  d.reps$ManEnvGen <- as.factor(paste(d.reps$Man, d.reps$Env, d.reps$Gen, sep = ":"))
  d.corr$ManEnvGen <- as.factor(paste(d.corr$Man, d.corr$Env, d.corr$Gen, sep = ":"))
  d.corr$reps <- numeric(nrow(d.corr))
  for (i in 1:nrow(d.corr)) {
    MEG <- as.character(d.corr$ManEnvGen)[i]
    reps <- length(droplevels(d.reps[d.reps$ManEnvGen == MEG,])$GY)
    d.corr[d.corr$ManEnvGen == MEG, "reps"] <- reps
  }
  
  # Loading in the kinship and subsetting to genotypes with SNPs:
  K <- readRDS("Briwecs/raw_data/K.rds")
  
  # Some genotypes have annoying names so we maually change those:
  rownames(K)[which(rownames(K) == "Triple dirk \"S\"")] <- "Triple Dirk S"
  rownames(K)[which(rownames(K) == "Boregan")] <- "Boregar"
  rownames(K)[which(rownames(K) == "Brillant")] <- "Brilliant"
  rownames(K)[which(rownames(K) == "Capelle Desprez")] <- "Cappelle_Desprez"
  rownames(K)[which(rownames(K) == "Lambriego Inia")] <- "Labriego-Inia"
  rownames(K)[which(rownames(K) == "Mex. 3")] <- "Mexico_3"
  rownames(K)[which(rownames(K) == "Robigous")] <- "Robigus"
  rownames(K)[which(rownames(K) == "Götz")] <- "Gotz"
  rownames(K)[which(rownames(K) == "Türkis")] <- "Turkis"
  rownames(K)[which(rownames(K) == "WW 4180 (Kongo)")] <- "WW 4180"
  
  colnames(K)[which(colnames(K) == "Triple dirk \"S\"")] <- "Triple Dirk S"
  colnames(K)[which(colnames(K) == "Boregan")] <- "Boregar"
  colnames(K)[which(colnames(K) == "Brillant")] <- "Brilliant"
  colnames(K)[which(colnames(K) == "Capelle Desprez")] <- "Cappelle_Desprez"
  colnames(K)[which(colnames(K) == "Lambriego Inia")] <- "Labriego-Inia"
  colnames(K)[which(colnames(K) == "Mex. 3")] <- "Mexico_3"
  colnames(K)[which(colnames(K) == "Robigous")] <- "Robigus"
  colnames(K)[which(colnames(K) == "Götz")] <- "Gotz"
  colnames(K)[which(colnames(K) == "Türkis")] <- "Turkis"
  colnames(K)[which(colnames(K) == "WW 4180 (Kongo)")] <- "WW 4180"
  
  rownames(K) <- gsub(" ", "_", rownames(K))
  colnames(K) <- gsub(" ", "_", colnames(K))
  
  d.corr <- droplevels(d.corr[d.corr$Gen %in% rownames(K),])
  d.corr$Man <- as.factor(d.corr$Man)
  d.corr$Env <- as.factor(d.corr$Env)
  d.corr$Gen <- as.factor(d.corr$Gen)
  K <- K[levels(d.corr$Gen), levels(d.corr$Gen)]
  saveRDS(d.corr, "Briwecs/data/d.corr.rds")
  saveRDS(K, "Briwecs/data/K.corr.rds")
} else {
  d.corr <- readRDS("Briwecs/data/d.corr.rds")
  K <- readRDS("Briwecs/data/K.corr.rds")
}

# For now, let's focus on LN/HN with no fungicide:
d.corr.ss <- droplevels(d.corr[d.corr$Man %in% c("LN_NF_RF", "HN_NF_RF"),])
d.corr.ss$Man <- as.factor(ifelse(as.character(d.corr.ss$Man) == "LN_NF_RF", "LN", "HN"))
d.corr.ss$ManEnv <- as.factor(paste(d.corr.ss$Man, d.corr.ss$Env, sep = ":"))

# Finding the common set of genotypes
tested <- vector("list", length(levels(d.corr.ss$ManEnv)))
for (i in 1:length(levels(d.corr.ss$ManEnv))) {
  tested[[i]] <- unique(as.character(droplevels(d.corr.ss[d.corr.ss$ManEnv == levels(d.corr.ss$ManEnv)[i],])$Gen))
  names(tested)[i] <- levels(d.corr.ss$ManEnv)[i]
}
tested <- Reduce(intersect, tested)
d.corr.ss <- droplevels(d.corr.ss[d.corr.ss$Gen %in% tested,])
K <- K[levels(d.corr.ss$Gen), levels(d.corr.ss$Gen)]

# ggplot(d.corr.ss, aes(x = Env, y = GY, color = Man)) + geom_boxplot()
# plot(droplevels(d.corr.ss[d.corr.ss$Man == "HN" & d.corr.ss$Env == "HAN_2017",])$GY,
#      droplevels(d.corr.ss[d.corr.ss$Man == "LN" & d.corr.ss$Env == "HAN_2017",])$GY)

fit <- asreml(fixed = GY ~ -1 + Env,
              random = ~ fa(Env, 3):vm(Gen, K),
              residual = ~ units,
              data = droplevels(d.corr.ss[d.corr.ss$Man == "HN",]))

PSI <- diag(summary(fit)$varcomp[1:14, "component"])
L <- matrix(0, 14, 3)
rownames(L) <- gsub(".*!(.*)!var", "\\1", rownames(summary(fit)$varcomp)[1:14])
L[, 1] <- summary(fit)$varcomp[15:28, "component"]
L[, 2] <- summary(fit)$varcomp[29:42, "component"]
L[, 3] <- summary(fit)$varcomp[43:56, "component"]

cormat.HN <- cov2cor(L %*% t(L) + PSI)
covmat.HN <- L %*% t(L) + PSI

fit <- asreml(fixed = GY ~ -1 + Env,
              random = ~ fa(Env, 3):vm(Gen, K),
              residual = ~ units,
              data = droplevels(d.corr.ss[d.corr.ss$Man == "LN",]))

PSI <- diag(summary(fit)$varcomp[1:14, "component"])
L <- matrix(0, 14, 3)
rownames(L) <- gsub(".*!(.*)!var", "\\1", rownames(summary(fit)$varcomp)[1:14])
L[, 1] <- summary(fit)$varcomp[15:28, "component"]
L[, 2] <- summary(fit)$varcomp[29:42, "component"]
L[, 3] <- summary(fit)$varcomp[43:56, "component"]

cormat.LN <- cov2cor(as.matrix(Matrix::nearPD(L %*% t(L) + PSI)$mat) + diag(10, 14, 14))
covmat.LN <- as.matrix(Matrix::nearPD(L %*% t(L) + PSI)$mat) + diag(10, 14, 14)

# Formatting and saving the benchmark data:
d.corr.ss$ManEnv <- as.factor(paste(d.corr.ss$Man, d.corr.ss$Env, sep = ":"))
d.corr.ss$EnvGen <- as.factor(paste(d.corr.ss$Env, d.corr.ss$Gen, sep = ":"))
d.corr.ss$ManGen <- as.factor(paste(d.corr.ss$Man, d.corr.ss$Gen, sep = ":"))
d.corr.ss$ManEnvGen <- as.factor(paste(d.corr.ss$Man, d.corr.ss$Env, d.corr.ss$Gen, sep = ":"))

M.man <- matrix(0, 2, 2); rownames(M.man) <- colnames(M.man) <- levels(d.corr.ss$Man)
M.env <- matrix(0, 14, 14); rownames(M.env) <- colnames(M.env) <- levels(d.corr.ss$Env)
K <- K[levels(d.corr.ss$Gen), levels(d.corr.ss$Gen)]

order <- rownames(kronecker(kronecker(M.man, M.env, make.dimnames = TRUE), K, make.dimnames = TRUE))

d.corr.ss <- d.corr.ss[match(order, as.character(d.corr.ss$ManEnvGen)),
                       c("Man", "Env", "Gen", "ManEnv", "ManGen", "EnvGen", "ManEnvGen", "GY", "reps")]

cormat.HN <- cormat.HN[levels(d.corr.ss$Env), levels(d.corr.ss$Env)]
cormat.LN <- cormat.LN[levels(d.corr.ss$Env), levels(d.corr.ss$Env)]

covmat.HN <- covmat.HN[levels(d.corr.ss$Env), levels(d.corr.ss$Env)]
covmat.LN <- covmat.LN[levels(d.corr.ss$Env), levels(d.corr.ss$Env)]

datalist.benchmark <- list(ydata = d.corr.ss,
                           cormat.HN = cormat.HN, cormat.LN = cormat.LN,
                           vars.HN = diag(covmat.HN), vars.LN = diag(covmat.LN))

saveRDS(datalist.benchmark, "Briwecs/data/datalist.benchmark.rds")
saveRDS(K, "Briwecs/data/K.benchmark.rds")

# Environmental data:
coords <- data.frame(Loc = c("HAN",      "KAL",      "KIE",      "QLB",     "RHH"),
                     Lat = c(52.2438183, 50.6131046, 54.3156846, 51.769213, 50.760497),
                     Lon = c(9.81810808, 6.9942357,  9.98041391, 11.145669, 8.875599))

d <- read.csv("Briwecs/raw_data/BRIWECS_data_publication.csv", sep = ";")
d$Env <- paste(d$Location, d$Year, sep = "_")
d <- droplevels(d[d$Env %in% rownames(datalist.benchmark$cormat.HN),])
d <- droplevels(d[d$Treatment %in% c("HN_NF_RF", "LN_NF_RF"),])
d$Treatment <- as.factor(ifelse(d$Treatment == "HN_NF_RF", "HN", "LN"))
d$ManEnv <- as.factor(paste(d$Treatment, d$Env, sep = ":"))
dates.sowing <- aggregate(d, Sowing_date ~ Env, FUN = function(x) round(mean(x)))
dates.harvest <- aggregate(d, BBCH87 ~ Env, FUN = function(x) round(mean(x)))

# Harvest dates are missing for KIE_2015, KIE_2016, RHH_2016, and RHH_2015.
# Best option I guess is to take the average trial duration, and add it to
# the sowing dates for these environments. Note that we're also averaging over
# managements because we want one set of env. covs. per environment.
dates.harvest <- rbind(dates.harvest,
                       expand.grid(Env = c("KIE_2015", "KIE_2016", "RHH_2015", "RHH_2016"),
                                   BBCH87 = NA))

dates.harvest <- dates.harvest[match(dates.sowing$Env, dates.harvest$Env),]
duration <- mean(na.omit(dates.harvest$BBCH87 + 365 - dates.sowing$Sowing_date))
dates.harvest[is.na(dates.harvest$BBCH87), "BBCH87"] <-
  round(dates.sowing[is.na(dates.harvest$BBCH87), "Sowing_date"] + duration - 365)

dates <- data.frame(Env = dates.sowing$Env)
dates$sowing <- as.Date(sprintf("%s-01-01", as.numeric(gsub(".*_(....)", "\\1", dates.sowing$Env)) - 1)) + dates.sowing$Sowing_date - 1
dates$harvest <- as.Date(sprintf("%s-01-01", as.numeric(gsub(".*_(....)", "\\1", dates.harvest$Env)))) + dates.harvest$BBCH87 - 1

dates$lon <- dates$lat <- rep(0, nrow(dates))
for (i in 1:nrow(coords)) {
  city <- coords$Loc[i]
  dates[grepl(city, dates$Env), c("lat", "lon")] <- coords[i, c("Lat", "Lon")]
}

df.clim <- EnvRtype::get_weather(env.id = dates$Env,
                                 lat = dates$lat,
                                 lon = dates$lon,
                                 start.day = dates$sowing,
                                 end.day = dates$harvest)

keep <- c("env",
          "DOY", "YYYYMMDD", "daysFromStart",
          "RH2M",                      # Relative humidity at 2m
          "WS2M",                      # Wind speed at 2m
          "PRECTOT",                   # Total precipitation
          "EVPTRNS", "P_ETP",          # Evapotranspiration
          "VPD",                       # Vapor pressure deficit
          "ALLSKY_SFC_LW_DWN",         # Radiation variables
          "ALLSKY_SFC_SW_DWN",
          "ALLSKY_SFC_SW_DNI",
          "ALLSKY_SFC_PAR_TOT",
          "ALLSKY_SFC_UVA",
          "ALLSKY_SFC_UVB",
          "RTA",                       # Extraterrestrial radiation
          "n", "N",                    # Actual duration of sunshine and daylight hours
          "T2M", "T2M_MAX", "T2M_MIN") # Mean, max, min temperature at 2m

edata <- df.clim[, keep]

# Checking whether we have data for all days from start to end:
for (i in unique(edata$env)) {
  days <- edata[edata$env == i, "daysFromStart"]
  if (all(days == seq(1, length(days)))) {
    cat(sprintf("%s is correct!\n", i))
  } else {
    cat(sprintf("%s is NOT correct!\n", i))
  }
}

# Now we calculate wheat growing degree days (https://ndawn.ndsu.nodak.edu/help-wheat-growing-degree-days.html):
edata$GDD <- edata$cumulGDD <- numeric(nrow(edata))
gdd <- function(tmin, tmax) {
  tmin <- tmin * 1.8 + 32 # Convert to F
  tmax <- tmax * 1.8 + 32
  if (tmin < 32) tmin <- 32 # Minimum and maximum growing temperatures for corn
  if (tmax < 32) tmax <- 32
  if (tmin > 95) tmin <- 95
  if (tmax > 95) tmax <- 95
  GDD <- ((tmin + tmax) / 2) - 32
  return(GDD)
}

# Calculate GDDs for individual dates:
for (i in 1:nrow(edata)) {
  edata$GDD[i] <- gdd(tmin = edata$T2M_MIN[i], tmax = edata$T2M_MAX[i])
}

# Calculate cumulative GDDs:
for (e in unique(edata$env)) {
  for (d in 1:nrow(edata[which(edata$env == e),])) {
    edata[which(edata$env == e & edata$day == d), "cumulGDD"] <-
      sum(edata[which(edata$env == e & edata$day %in% 1:d), "GDD"])
  }
}

# Printing end of season GDDs (not that much variation which is good):
for (e in unique(edata$env)) {
  cat(sprintf("Environment: %s\tEOF GDD: %s\n", e, round(max(edata[which(edata$env == e), "cumulGDD"]))))
}

vars <- colnames(edata)[5:22]

# Cumulative GDD windows:
wins <- c("(0,100]", "(100,200]", "(200,300]", "(300,400]", "(400,500]", "(500,600]", "(600,700]", "(700,800]",
          "(800,900]", "(900,1e+03]", "(1e+03,1.1e+03]", "(1.1e+03,1.2e+03]", "(1.2e+03,1.3e+03]",
          "(1.3e+03,1.4e+03]", "(1.4e+03,1.5e+03]", "(1.5e+03,1.6e+03]", "(1.6e+03,1.7e+03]", "(1.7e+03,1.8e+03]",
          "(1.8e+03,1.9e+03]", "(1.9e+03,2e+03]", "(2e+03,2.1e+03]", "(2.1e+03,2.2e+03]", "(2.2e+03,2.3e+03]",
          "(2.3e+03,2.4e+03]", "(2.4e+03,2.5e+03]", "(2.5e+03,2.6e+03]", "(2.6e+03,2.7e+03]", "(2.7e+03,2.8e+03]",
          "(2.8e+03,2.9e+03]", "(2.9e+03,3e+03]", "(3e+03,3.1e+03]", "(3.1e+03,3.2e+03]", "(3.2e+03,3.3e+03]",
          "(3.3e+03,3.4e+03]", "(3.4e+03,3.5e+03]", "(3.5e+03,3.6e+03]", "(3.6e+03,3.7e+03]", "(3.7e+03,3.8e+03]",
          "(3.8e+03,3.9e+03]", "(3.9e+03,4.0e+03]", "(4.0e+03,4.1e+03]", "(4.1e+03,4.2e+03]", "(4.2e+03,4.3e+03]")

edata.fit <- expand.grid(env = unique(edata$env),
                         var = vars,
                         window = wins,
                         AUC = NA,
                         stringsAsFactors = F)

# Scaling environmental variables:
edata[, 5:22] <- scale(edata[, 5:22])

# Fitting piecewiese constant regressions for each variable with breakpoints at each
# transition from one window to the next.
for (e in unique(edata.fit$env)) {
  for (t in unique(edata.fit$var)) {
    
    edata.ss <- na.omit(edata[which(edata$env == e), c(t, "cumulGDD")])
    edata.ss$window <- cut(edata.ss$cumulGDD, breaks = seq(0, 4300, 100))
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
datalist <- list(ydata = d.corr.ss)
X <- scale(t(as.matrix(edata.fit.wide[, -1])))
C <- (t(X) %*% X) / (ncol(edata.fit.wide[, -1]) - 1)
C <- C[rownames(datalist.benchmark$cormat.HN), rownames(datalist.benchmark$cormat.HN)]
datalist$EC <- C

# Distance matrix:
X <- scale(as.matrix(edata.fit.wide[, -1]))
ED <- (as.matrix(dist(X, method = "euclidian"))^2) / ncol(X)
ED <- ED[rownames(datalist.benchmark$cormat.HN), colnames(datalist.benchmark$cormat.HN)]
datalist$ED <- ED

# Trying some stuff:
grid <- seq(0.01, 2, length.out = 500)
diff <- function(R, D, h) {
  Rh <- exp(-h * D)
  return(mean(abs(R - Rh)))
}

d <- numeric(length(grid))
cormat <- (datalist.benchmark$cormat.LN + datalist.benchmark$cormat.HN) / 2
for (h in 1:length(grid)) {
  d[h] <- diff(R = cormat, D = ED, h = grid[h])
}
plot(y = d, x = grid)
(h.opt <- grid[which(d == min(d))])
round(cormat - exp(-h.opt * ED), 2)
mean(abs(cormat[upper.tri(cormat)] - exp(-h.opt * ED)[upper.tri(cormat)]))

all(unique(datalist$ydata$Man) == levels(datalist$ydata$Man))
all(unique(datalist$ydata$Env) == levels(datalist$ydata$Env))
all(unique(datalist$ydata$Gen) == levels(datalist$ydata$Gen))

all(levels(datalist$ydata$Env) == rownames(datalist$EC))
all(levels(datalist$ydata$Gen) == rownames(K))
all(levels(datalist$ydata$EnvGen) == rownames(kronecker(datalist$EC, K, make.dimnames = TRUE)))

saveRDS(datalist, "Briwecs/data/datalist.rds")
saveRDS(K, "Briwecs/data/K.rds")

d.full <- droplevels(d.full[d.full$Treatment %in% c("LN_NF_RF", "HN_NF_RF"),])
d.full$Man <- as.factor(ifelse(as.character(d.full$Treatment) == "LN_NF_RF", "LN", "HN"))
d.full$ManEnvGen <- as.factor(paste(d.full$Man, d.full$Env, d.full$BRISONr, sep = ":"))
d.misc <- aggregate(d.full, cbind(Sowing_date, Emergence_date, BBCH59, BBCH87,
                                  Plantheight_bio, Seedyield, Seedyield_bio, Biomass_bio,
                                  Harvest_Index_bio, TGW, TGW_bio, Spike_number_bio,
                                  Stripe_rust, Powdery_mildew, Leaf_rust, Septoria,
                                  DTR, Fusarium, Falling_number, Crude_protein,
                                  Sedimentation, Grain_per_spike_bio, Grain, Biomass,
                                  Protein_yield
                                  ) ~ ManEnvGen, FUN = mean, na.action = na.pass)
d.misc <- droplevels(d.misc[match(datalist$ydata$ManEnvGen, d.misc$ManEnvGen),])

info <- read.csv("Briwecs/raw_data/BRIWECs_cultivar_info.csv")
# Some genotypes have annoying names so we maually change those:
info$genotype[which(info$genotype == "Triple Dirk \"S\"")] <- "Triple Dirk S"
info$genotype[which(info$genotype == "G\xf6tz")] <- "Gotz"
info$genotype[which(info$genotype == "T\xfcrkis")] <- "Turkis"

d.misc$Gen <- as.factor(gsub("..:.*:(.*)", "\\1", as.character(d.misc$ManEnvGen)))
info$genotype <- gsub(" ", "_", info$genotype)
info <- droplevels(info[match(d.misc$Gen, info$genotype),])

d.misc <- cbind(d.misc, info[, c("baking_qulaity", "country", "breeder", "RYear")])
datalist$ydata <- cbind(datalist$ydata, d.misc[, setdiff(colnames(d.misc), c("ManEnvGen", "Gen"))])

saveRDS(datalist, "Briwecs/data/datalist.miscinfo.rds")

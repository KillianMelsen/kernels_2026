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

# Loading and subsetting yield data:
ydata <- read.csv("DROPS/raw_data/2b-GrainYield_components_BLUEs_level-1.csv")
ydata <- ydata[which(ydata$Variety_ID %in% rownames(K)),]

# Find the trials that have all 246 varieties and subset to those:
keep <- names(table(ydata$Experiment))[which(table(ydata$Experiment) == 246)]
ydata <- ydata[which(ydata$Experiment %in% keep),]; rm(keep)
ydata <- droplevels(ydata[ydata$Experiment != "Deb13R",]) # Does not make a lot of difference to exclude it for kernel computation already.

# Loading environmental data:
edata <- read.csv("DROPS/raw_data/1-Env_variables_daily.csv")

# Which loc-years do we have in ydata (14 loc-years)?
y.locyears <- unique(gsub("(.*)[RW]", "\\1", unique(ydata$Experiment)))

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
saveRDS(C, "DROPS/data/EC.rds")

# Distance matrix:
X <- scale(as.matrix(edata.fit.wide[, -1]))
ED <- (as.matrix(dist(X, method = "euclidian"))^2) / ncol(X)
saveRDS(ED, "DROPS/data/ED.rds")

# Subsetting data to rainfed trials:
keep <- c(paste0(y.locyears, "R"), paste0(y.locyears, "W"))
ydata.ss <- droplevels(ydata[which(ydata$Experiment %in% keep),])
ydata.ss$Env <- substr(ydata.ss$Experiment, 1, 5)
ydata.ss$Man <- substr(ydata.ss$Experiment, 6, 6)
ydata.ss <- ydata.ss[, c("Env", "Man", "Variety_ID", "grain.yield", "grain.number",
                         "seed.size", "plant.height", "tassel.height", "ear.height",
                         "anthesis", "silking")]
names(ydata.ss) <- c("Env", "Man", "Gen", "GY", "GN", "SS", "PH", "TH", "EH", "A", "S")

# Making 100% sure the ordering is right:
M <- matrix(0, 2, 2)
rownames(M) <- colnames(M) <- c("R", "W")
order <- rownames(kronecker(kronecker(M, ED, make.dimnames = TRUE), K, make.dimnames = TRUE))
ydata.ss$ManEnvGen <- paste(ydata.ss$Man, ydata.ss$Env, ydata.ss$Gen, sep = ":")
ydata.ss$ManEnv <- paste(ydata.ss$Man, ydata.ss$Env, sep = ":")
ydata.ss$ManGen <- paste(ydata.ss$Man, ydata.ss$Gen, sep = ":")
ydata.ss$EnvGen <- paste(ydata.ss$Env, ydata.ss$Gen, sep = ":")
ydata.ss <- ydata.ss[match(order, ydata.ss$ManEnvGen),]

ydata.ss$Env <- factor(ydata.ss$Env, levels = ydata.ss$Env, labels = ydata.ss$Env)
ydata.ss$Man <- factor(ydata.ss$Man, levels = ydata.ss$Man, labels = ydata.ss$Man)
ydata.ss$Gen <- factor(ydata.ss$Gen, levels = ydata.ss$Gen, labels = ydata.ss$Gen)
ydata.ss$ManEnvGen <- factor(ydata.ss$ManEnvGen, levels = ydata.ss$ManEnvGen, labels = ydata.ss$ManEnvGen)
ydata.ss$ManEnv <- factor(ydata.ss$ManEnv, levels = ydata.ss$ManEnv, labels = ydata.ss$ManEnv)
ydata.ss$ManGen <- factor(ydata.ss$ManGen, levels = ydata.ss$ManGen, labels = ydata.ss$ManGen)
ydata.ss$EnvGen <- factor(ydata.ss$EnvGen, levels = ydata.ss$EnvGen, labels = ydata.ss$EnvGen)

all(levels(ydata.ss$Env) == rownames(ED))
all(levels(ydata.ss$Gen) == rownames(K))
all(levels(ydata.ss$ManEnv) == rownames(kronecker(M, ED, make.dimnames = TRUE)))
all(levels(ydata.ss$ManEnvGen) == rownames(kronecker(kronecker(M, ED, make.dimnames = TRUE), K, make.dimnames = TRUE)))

all(levels(ydata.ss$Man) == unique(as.character(ydata.ss$Man)))
all(levels(ydata.ss$Env) == unique(as.character(ydata.ss$Env)))
all(levels(ydata.ss$Gen) == unique(as.character(ydata.ss$Gen)))

all(levels(ydata.ss$ManEnvGen) == unique(as.character(ydata.ss$ManEnvGen)))

all(levels(ydata.ss$ManEnv) == unique(as.character(ydata.ss$ManEnv)))
all(levels(ydata.ss$ManGen) == unique(as.character(ydata.ss$ManGen)))
all(levels(ydata.ss$EnvGen) == unique(as.character(ydata.ss$EnvGen)))

# Saving:
saveRDS(ydata.ss, "DROPS/data/ydata.rds")



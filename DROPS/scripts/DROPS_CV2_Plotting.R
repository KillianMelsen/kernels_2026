#!/usr/bin/env Rscript

library(ggplot2)
library(patchwork)
results <- readRDS("DROPS/results/CV2/results.CV2.SE.1-150.rds")
results.mb <- readRDS("DROPS/results/CV2/old_MB_results/results.MB.CV2.SE.1-150.rds")
results <- rbind(results, results.mb)
means <- aggregate(results, cbind(cor_pearson, RMSE) ~ Checks + Model + Man, FUN = mean)
SEs <- aggregate(results, cbind(cor_pearson, RMSE) ~ Checks + Model + Man, FUN = function(x) sd(x)/sqrt(length(x)))
names(SEs)[4:5] <- c("SE_cor_pearson", "SE_RMSE")
plotdata <- cbind(means, SEs[, 4:5])
plotdata.corr <- plotdata[, c(1:3, 4, 6)]
plotdata.rmse <- plotdata[, c(1:3, 5, 7)]
plotdata.corr$Measure <- "Pearson Correlation"
plotdata.rmse$Measure <- "RMSE"
names(plotdata.corr)[4] <- names(plotdata.rmse)[4] <- "Value"
names(plotdata.corr)[5] <- names(plotdata.rmse)[5] <- "SE"
plotdata.long <- rbind(plotdata.corr, plotdata.rmse)

plotdata.long$Model <- factor(plotdata.long$Model,
                              levels = levels(plotdata.long$Model),
                              labels = c("ADD", "FA-1", "FA-2", "FA-3",
                                         "SV-LK", "SV-GK", "MV-LK", "MV-GK", "MB-SV-GK", "MB-MV-GK"))
plotdata.long$Man <- factor(plotdata.long$Man, levels = levels(plotdata.long$Man), labels = c("Rain-fed", "Irrigated"))

# Text of section 3.1 ====
bw <- readRDS("DROPS/results/CV2/bandwidths.MB.CV2.SE.1-150.rds")
aggregate(bw, cbind(R, W, R_W) ~ Model, FUN = function(x) round(mean(x), 2))

p1 <- ggplot(droplevels(plotdata.long[plotdata.long$Measure == "Pearson Correlation",]), aes(x = Checks, y = Value, color = Model)) +
  facet_wrap(vars(Man), ncol = 2) +
  geom_ribbon(aes(y = Value, ymin = Value - SE, ymax = Value + SE, fill = Model), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 2) +
  theme_classic(base_size = 18) +
  ylab("Pearson Correlation") + #ylim(c(0.545, 0.70)) +
  xlab(NULL) +
  scale_color_manual(values = c("black", "#0C62AF", "#4499F5", "#8FCAFD", "#FD8700", "#D8511D", "#758717", "#42673B", "pink1", "purple")) +
  scale_fill_manual(values = c("black", "#0C62AF", "#4499F5", "#8FCAFD", "#FD8700", "#D8511D", "#758717", "#42673B", "pink1", "purple")) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 20),
        axis.title.y = element_text(size = 20))

p2 <- ggplot(droplevels(plotdata.long[plotdata.long$Measure == "RMSE",]), aes(x = Checks, y = Value, color = Model)) +
  facet_wrap(vars(Man), ncol = 2, scales = "free") +
  geom_ribbon(aes(y = Value, ymin = Value - SE, ymax = Value + SE, fill = Model), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 2) +
  theme_classic(base_size = 18) +
  ylab("RMSE") +
  xlab("Number of checks") +
  scale_color_manual(values = c("black", "#0C62AF", "#4499F5", "#8FCAFD", "#FD8700", "#D8511D", "#758717", "#42673B", "pink1", "purple")) +
  scale_fill_manual(values = c("black", "#0C62AF", "#4499F5", "#8FCAFD", "#FD8700", "#D8511D", "#758717", "#42673B", "pink1", "purple")) +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y = element_text(size = 20))

p1 / p2
ggsave(filename = "plots/DROPS_CV2.png", dpi = 300, width = 32, height = 32, units = "cm")
# rm(list = ls())

# Enlarged versions:
PCR <- ggplot(droplevels(plotdata.long[plotdata.long$Measure == "Pearson Correlation" & plotdata.long$Man == "Rain-fed",]), aes(x = Checks, y = Value, color = Model)) +
  facet_wrap(vars(Man), ncol = 2) +
  geom_ribbon(aes(y = Value, ymin = Value - SE, ymax = Value + SE, fill = Model), alpha = 0.2, color = NA) +
  geom_line(linewidth = 2.5) +
  geom_point(size = 3) +
  theme_classic(base_size = 18) +
  ylab("Pearson Correlation") + #ylim(c(0.545, 0.70)) +
  xlab("Number of checks") +
  scale_color_manual(values = c("black", "#0C62AF", "#4499F5", "#8FCAFD", "#FD8700", "#D8511D", "#758717", "#42673B", "pink1", "purple")) +
  scale_fill_manual(values = c("black", "#0C62AF", "#4499F5", "#8FCAFD", "#FD8700", "#D8511D", "#758717", "#42673B", "pink1", "purple")) +
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 20),
        axis.title.y = element_text(size = 20))
ggsave(plot = PCR, filename = "plots/DROPS_CV2_PCR.png", dpi = 300, width = 32, height = 48, units = "cm")

PCW <- ggplot(droplevels(plotdata.long[plotdata.long$Measure == "Pearson Correlation" & plotdata.long$Man == "Irrigated",]), aes(x = Checks, y = Value, color = Model)) +
  facet_wrap(vars(Man), ncol = 2) +
  geom_ribbon(aes(y = Value, ymin = Value - SE, ymax = Value + SE, fill = Model), alpha = 0.2, color = NA) +
  geom_line(linewidth = 2.5) +
  geom_point(size = 3) +
  theme_classic(base_size = 18) +
  ylab("Pearson Correlation") + #ylim(c(0.545, 0.70)) +
  xlab("Number of checks") +
  scale_color_manual(values = c("black", "#0C62AF", "#4499F5", "#8FCAFD", "#FD8700", "#D8511D", "#758717", "#42673B", "pink1", "purple")) +
  scale_fill_manual(values = c("black", "#0C62AF", "#4499F5", "#8FCAFD", "#FD8700", "#D8511D", "#758717", "#42673B", "pink1", "purple")) +
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 20),
        axis.title.y = element_text(size = 20))
ggsave(plot = PCW, filename = "plots/DROPS_CV2_PCW.png", dpi = 300, width = 32, height = 48, units = "cm")

RMSER <- ggplot(droplevels(plotdata.long[plotdata.long$Measure == "RMSE" & plotdata.long$Man == "Rain-fed",]), aes(x = Checks, y = Value, color = Model)) +
  facet_wrap(vars(Man), ncol = 2) +
  geom_ribbon(aes(y = Value, ymin = Value - SE, ymax = Value + SE, fill = Model), alpha = 0.2, color = NA) +
  geom_line(linewidth = 2.5) +
  geom_point(size = 3) +
  theme_classic(base_size = 18) +
  ylab("RMSE") + #ylim(c(0.545, 0.70)) +
  xlab("Number of checks") +
  scale_color_manual(values = c("black", "#0C62AF", "#4499F5", "#8FCAFD", "#FD8700", "#D8511D", "#758717", "#42673B", "pink1", "purple")) +
  scale_fill_manual(values = c("black", "#0C62AF", "#4499F5", "#8FCAFD", "#FD8700", "#D8511D", "#758717", "#42673B", "pink1", "purple")) +
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 20),
        axis.title.y = element_text(size = 20))
ggsave(plot = RMSER, filename = "plots/DROPS_CV2_RMSER.png", dpi = 300, width = 32, height = 48, units = "cm")

RMSEW <- ggplot(droplevels(plotdata.long[plotdata.long$Measure == "RMSE" & plotdata.long$Man == "Irrigated",]), aes(x = Checks, y = Value, color = Model)) +
  facet_wrap(vars(Man), ncol = 2) +
  geom_ribbon(aes(y = Value, ymin = Value - SE, ymax = Value + SE, fill = Model), alpha = 0.2, color = NA) +
  geom_line(linewidth = 2.5) +
  geom_point(size = 3) +
  theme_classic(base_size = 18) +
  ylab("RMSE") + #ylim(c(0.545, 0.70)) +
  xlab("Number of checks") +
  scale_color_manual(values = c("black", "#0C62AF", "#4499F5", "#8FCAFD", "#FD8700", "#D8511D", "#758717", "#42673B", "pink1", "purple")) +
  scale_fill_manual(values = c("black", "#0C62AF", "#4499F5", "#8FCAFD", "#FD8700", "#D8511D", "#758717", "#42673B", "pink1", "purple")) +
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 20),
        axis.title.y = element_text(size = 20))
ggsave(plot = RMSEW, filename = "plots/DROPS_CV2_RMSEW.png", dpi = 300, width = 32, height = 48, units = "cm")


# Plotting per environment:
rm(list = ls())
results <- readRDS("DROPS/results/CV2/results.CV2.SE.1-150.rds")
results.mb <- readRDS("DROPS/results/CV2/old_MB_results/results.MB.CV2.SE.1-150.rds")
results <- rbind(results, results.mb)
means <- aggregate(results, cbind(cor_pearson, RMSE) ~ Checks + Model + Man + Env, FUN = mean)
SEs <- aggregate(results, cbind(cor_pearson, RMSE) ~ Checks + Model + Man + Env, FUN = function(x) sd(x)/sqrt(length(x)))
names(SEs)[5:6] <- c("SE_cor_pearson", "SE_RMSE")
plotdata <- cbind(means, SEs[, 5:6])
plotdata.corr <- plotdata[, c(1:4, 5, 7)]
plotdata.rmse <- plotdata[, c(1:4, 6, 8)]
plotdata.corr$Measure <- "Pearson Correlation"
plotdata.rmse$Measure <- "RMSE"
names(plotdata.corr)[5] <- names(plotdata.rmse)[5] <- "Value"
names(plotdata.corr)[6] <- names(plotdata.rmse)[6] <- "SE"
plotdata.long <- rbind(plotdata.corr, plotdata.rmse)

plotdata.long$Model <- factor(plotdata.long$Model,
                              levels = levels(plotdata.long$Model),
                              labels = c("ADD", "FA-1", "FA-2", "FA-3",
                                         "SV-LK", "SV-GK", "MV-LK", "MV-GK", "MB-SV-GK", "MB-MV-GK"))
plotdata.long$Man <- factor(plotdata.long$Man, levels = levels(plotdata.long$Man), labels = c("Rain-fed", "Irrigated"))
# plotdata.long.corr.R <- droplevels(plotdata.long[plotdata.long$Measure == "Pearson correlation" & plotdata.long$Man == "Rain-fed",])
# plotdata.long.corr.W <- droplevels(plotdata.long[plotdata.long$Measure == "Pearson correlation" & plotdata.long$Man == "Irrigated",])
# plotdata.long.rmse.R <- droplevels(plotdata.long[plotdata.long$Measure == "RMSE" & plotdata.long$Man == "Rain-fed",])
# plotdata.long.rmse.W <- droplevels(plotdata.long[plotdata.long$Measure == "RMSE" & plotdata.long$Man == "Irrigated",])

for (i in 1:4) {
  if (i %in% c(1, 3)) {
    envs <- levels(plotdata.long$Env)[1:7]
  } else {
    envs <- levels(plotdata.long$Env)[8:14]
  }
  metr <- ifelse(i %in% c(1, 2), "Pearson Correlation", "RMSE")
  ggplot(droplevels(plotdata.long[plotdata.long$Env %in% envs & plotdata.long$Measure == metr,]), aes(x = Checks, y = Value, color = Model)) +
    facet_grid(rows = vars(Env), cols = vars(Man), scales = "free_y") +
    geom_ribbon(aes(y = Value, ymin = Value - SE, ymax = Value + SE, fill = Model), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.5) +
    geom_point(size = 2) +
    theme_classic(base_size = 20) +
    ylab(metr) +
    xlab("Number of checks") +
    scale_color_manual(values = c("black", "#0C62AF", "#4499F5", "#8FCAFD", "#FD8700", "#D8511D", "#758717", "#42673B", "pink1", "purple")) +
    scale_fill_manual(values = c("black", "#0C62AF", "#4499F5", "#8FCAFD", "#FD8700", "#D8511D", "#758717", "#42673B", "pink1", "purple")) +
    theme(legend.position = "bottom",
          strip.background = element_blank(),
          strip.text = element_text(size = 20),
          axis.title.y = element_text(size = 20))
  
  label1 <- ifelse(i %in% c(1, 2), "CORR", "RMSE")
  label2 <- ifelse(i %in% c(1, 3), "A", "B")
  ggsave(filename = sprintf("plots/accuracies_per_env/DROPS_%s_%s.png", label1, label2), dpi = 300, width = 32, height = 48, units = "cm")
}



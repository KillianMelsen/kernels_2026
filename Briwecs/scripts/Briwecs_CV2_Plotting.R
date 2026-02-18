#!/usr/bin/env Rscript

library(ggplot2)
library(patchwork)
results <- readRDS("Briwecs/results/CV2/results.CV2.weights.1-250.rds")
means <- aggregate(results, cbind(cor_pearson, RMSE) ~ Checks + Model + Man, FUN = mean)
SEs <- aggregate(results, cbind(cor_pearson, RMSE) ~ Checks + Model + Man, FUN = function(x) sd(x)/sqrt(length(x)))
names(SEs)[4:5] <- c("SE_cor_pearson", "SE_RMSE")
plotdata <- cbind(means, SEs[, 4:5])
plotdata.corr <- plotdata[, c(1:3, 4, 6)]
plotdata.rmse <- plotdata[, c(1:3, 5, 7)]
plotdata.corr$Measure <- "Pearson correlation"
plotdata.rmse$Measure <- "RMSE"
names(plotdata.corr)[4] <- names(plotdata.rmse)[4] <- "Value"
names(plotdata.corr)[5] <- names(plotdata.rmse)[5] <- "SE"
plotdata.long <- rbind(plotdata.corr, plotdata.rmse)

plotdata.long$Model <- factor(plotdata.long$Model,
                              levels = levels(plotdata.long$Model),
                              labels = c("ME", "FA-1", "FA-2", "FA-3",
                                         "SV-LK", "SV-GK", "MV-LK", "MV-GK"))
plotdata.long$Man <- factor(plotdata.long$Man, levels = levels(plotdata.long$Man), labels = c("High Nitrogen", "Low Nitrogen"))
ggplot(plotdata.long, aes(x = Checks, y = Value, color = Model)) +
  facet_wrap(vars(Measure, Man), ncol = 2, scales = "free") +
  geom_ribbon(aes(y = Value, ymin = Value - SE, ymax = Value + SE, fill = Model), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 2) +
  theme_classic(base_size = 18) +
  ylab(NULL) +
  xlab("Number of checks") +
  scale_color_manual(values = c("black", "#0C62AF", "#4499F5", "#8FCAFD", "#FD8700", "#D8511D", "#758717", "#42673B")) +
  scale_fill_manual(values = c("black", "#0C62AF", "#4499F5", "#8FCAFD", "#FD8700", "#D8511D", "#758717", "#42673B")) +
  theme(legend.position = "bottom")

ggsave(filename = "plots/Briwecs_CV2.png", dpi = 300, width = 32, height = 32, units = "cm")

p1 <- ggplot(droplevels(plotdata.long[plotdata.long$Measure == "Pearson correlation",]), aes(x = Checks, y = Value, color = Model)) +
  facet_wrap(vars(Measure, Man), ncol = 2) +
  geom_ribbon(aes(y = Value, ymin = Value - SE, ymax = Value + SE, fill = Model), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 2) +
  theme_classic(base_size = 18) +
  ylab(NULL) +
  xlab(NULL) +
  scale_color_manual(values = c("black", "#0C62AF", "#4499F5", "#8FCAFD", "#FD8700", "#D8511D", "#758717", "#42673B")) +
  scale_fill_manual(values = c("black", "#0C62AF", "#4499F5", "#8FCAFD", "#FD8700", "#D8511D", "#758717", "#42673B")) +
  theme(legend.position = "none")

p2 <- ggplot(droplevels(plotdata.long[plotdata.long$Measure == "RMSE",]), aes(x = Checks, y = Value, color = Model)) +
  facet_wrap(vars(Measure, Man), ncol = 2) +
  geom_ribbon(aes(y = Value, ymin = Value - SE, ymax = Value + SE, fill = Model), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 2) +
  theme_classic(base_size = 18) +
  ylab(NULL) +
  xlab("Number of checks") +
  scale_color_manual(values = c("black", "#0C62AF", "#4499F5", "#8FCAFD", "#FD8700", "#D8511D", "#758717", "#42673B")) +
  scale_fill_manual(values = c("black", "#0C62AF", "#4499F5", "#8FCAFD", "#FD8700", "#D8511D", "#758717", "#42673B")) +
  theme(legend.position = "bottom")

p1 / p2
ggsave(filename = "plots/Briwecs_CV2_V2.png", dpi = 300, width = 32, height = 32, units = "cm")

# rm(list = ls())


# Plotting per environment:
rm(list = ls())
results <- readRDS("Briwecs/results/CV2/results.CV2.1-250.rds")
means <- aggregate(results, cbind(cor_pearson, RMSE) ~ Checks + Model + Man + Env, FUN = mean)
SEs <- aggregate(results, cbind(cor_pearson, RMSE) ~ Checks + Model + Man + Env, FUN = function(x) sd(x)/sqrt(length(x)))
names(SEs)[5:6] <- c("SE_cor_pearson", "SE_RMSE")
plotdata <- cbind(means, SEs[, 5:6])
plotdata.corr <- plotdata[, c(1:4, 5, 7)]
plotdata.rmse <- plotdata[, c(1:4, 6, 8)]
plotdata.corr$Measure <- "Pearson correlation"
plotdata.rmse$Measure <- "RMSE"
names(plotdata.corr)[5] <- names(plotdata.rmse)[5] <- "Value"
names(plotdata.corr)[6] <- names(plotdata.rmse)[6] <- "SE"
plotdata.long <- rbind(plotdata.corr, plotdata.rmse)

plotdata.long$Model <- factor(plotdata.long$Model,
                              levels = levels(plotdata.long$Model),
                              labels = c("ME", "FA-1", "FA-2", "FA-3",
                                         "SV-LK", "SV-GK", "MV-LK", "MV-GK"))
plotdata.long$Man <- factor(plotdata.long$Man, levels = levels(plotdata.long$Man), labels = c("High Nitrogen", "Low Nitrogen"))
# plotdata.long.corr.HN <- droplevels(plotdata.long[plotdata.long$Measure == "Pearson correlation" & plotdata.long$Man == "High nitrogen",])
# plotdata.long.corr.LN <- droplevels(plotdata.long[plotdata.long$Measure == "Pearson correlation" & plotdata.long$Man == "Low nitrogen",])
# plotdata.long.rmse.HN <- droplevels(plotdata.long[plotdata.long$Measure == "RMSE" & plotdata.long$Man == "High nitrogen",])
# plotdata.long.rmse.LN <- droplevels(plotdata.long[plotdata.long$Measure == "RMSE" & plotdata.long$Man == "Low nitrogen",])

for (i in 1:4) {
  if (i %in% c(1, 3)) {
    envs <- levels(plotdata.long$Env)[1:7]
  } else {
    envs <- levels(plotdata.long$Env)[8:14]
  }
  metr <- ifelse(i %in% c(1, 2), "Pearson correlation", "RMSE")
  ggplot(droplevels(plotdata.long[plotdata.long$Env %in% envs & plotdata.long$Measure == metr,]), aes(x = Checks, y = Value, color = Model)) +
    facet_grid(rows = vars(Env), cols = vars(Man), scales = "free_y") +
    geom_ribbon(aes(y = Value, ymin = Value - SE, ymax = Value + SE, fill = Model), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.5) +
    geom_point(size = 2) +
    theme_classic(base_size = 20) +
    ylab(metr) +
    xlab("Number of checks") +
    scale_color_manual(values = c("black", "#0C62AF", "#4499F5", "#8FCAFD", "#FD8700", "#D8511D", "#758717", "#42673B")) +
    scale_fill_manual(values = c("black", "#0C62AF", "#4499F5", "#8FCAFD", "#FD8700", "#D8511D", "#758717", "#42673B")) +
    theme(legend.position = "bottom"#,
          # panel.background = element_rect(fill = 'transparent'), #transparent panel bg
          # plot.background = element_rect(fill='transparent', color = NA),
          # legend.background = element_rect(fill='transparent')
    )
  
  label1 <- ifelse(i %in% c(1, 2), "CORR", "RMSE")
  label2 <- ifelse(i %in% c(1, 3), "A", "B")
  ggsave(filename = sprintf("plots/accuracies_per_env/BRIWECS_%s_%s.png", label1, label2), dpi = 300, width = 32, height = 48, units = "cm")
}





























# plotdata.long.corr <- droplevels(plotdata.long[plotdata.long$Measure == "Pearson correlation",])
# plotdata.long.rmse <- droplevels(plotdata.long[plotdata.long$Measure == "RMSE",])
# 
# ggplot(plotdata.long.corr, aes(x = Checks, y = Value, color = Model)) +
#   facet_wrap(vars(Man), ncol = 2, scales = "free") +
#   geom_ribbon(aes(y = Value, ymin = Value - SE, ymax = Value + SE, fill = Model), alpha = 0.2, color = NA) +
#   geom_line(linewidth = 1.5) +
#   geom_point(size = 2) +
#   theme_classic(base_size = 20) +
#   ylab("Pearson correlation") +
#   xlab("Number of checks") +
#   scale_color_manual(values = c("black", "#0C62AF", "#4499F5", "#8FCAFD", "#FD8700", "#D8511D", "#758717", "#42673B")) +
#   scale_fill_manual(values = c("black", "#0C62AF", "#4499F5", "#8FCAFD", "#FD8700", "#D8511D", "#758717", "#42673B")) +
#   theme(legend.position = "bottom",
#         panel.background = element_rect(fill = 'transparent'), #transparent panel bg
#         plot.background = element_rect(fill='transparent', color=NA),
#         legend.background = element_rect(fill='transparent'))
# 
# ggsave(filename = "plots/Briwecs_CV2_corr.png", dpi = 300, width = 32, height = 20, units = "cm", bg = "transparent")
# 
# ggplot(plotdata.long.rmse, aes(x = Checks, y = Value, color = Model)) +
#   facet_wrap(vars(Man), ncol = 2, scales = "free") +
#   geom_ribbon(aes(y = Value, ymin = Value - SE, ymax = Value + SE, fill = Model), alpha = 0.2, color = NA) +
#   geom_line(linewidth = 1.5) +
#   geom_point(size = 2) +
#   theme_classic(base_size = 20) +
#   ylab("RMSE") +
#   xlab("Number of checks") +
#   scale_color_manual(values = c("black", "#0C62AF", "#4499F5", "#8FCAFD", "#FD8700", "#D8511D", "#758717", "#42673B")) +
#   scale_fill_manual(values = c("black", "#0C62AF", "#4499F5", "#8FCAFD", "#FD8700", "#D8511D", "#758717", "#42673B")) +
#   theme(legend.position = "bottom",
#         panel.background = element_rect(fill = 'transparent'), #transparent panel bg
#         plot.background = element_rect(fill='transparent', color=NA),
#         legend.background = element_rect(fill='transparent'))
# 
# ggsave(filename = "plots/Briwecs_CV2_rmse.png", dpi = 300, width = 32, height = 20, units = "cm", background = "transparent")




# Plotting convergence
# results <- readRDS("Briwecs/results/CV2/results.CV2.1-30.rds")
# means <- aggregate(results, Converged ~ Model + Checks, FUN = mean)
# 
# ggplot(means, aes(x = Checks, y = Converged * 100, color = Model)) +
#   geom_line(linewidth = 1.5) +
#   geom_point(size = 2) +
#   theme_classic() +
#   ylab("Percentage of converged fits (25 iterations)") +
#   xlab("Number of checks") +
#   scale_color_manual(values = c("black", "#0C62AF", "#4499F5", "#8FCAFD", "#FD8700", "#D8511D", "#758717", "#42673B")) +
#   scale_y_continuous(labels = scales::unit_format(unit = "%"))



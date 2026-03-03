# Plotting the locations ====
library(sf)
library(ggplot2)
library(rnaturalearth)
library(ggmap)

color = "black"

# Briwecs data
coords <- data.frame(Loc = c("HAN",      "KAL",      "KIE",      "QLB",     "RHH"),
                     Lat = c(52.2438183, 50.6131046, 54.3156846, 51.769213, 50.760497),
                     Lon = c(9.81810808, 6.9942357,  9.98041391, 11.145669, 8.875599))

DE_map <- rnaturalearth::ne_countries(country = "germany", returnclass = "sf")
# DE_map <-  get_stadiamap(bbox = unname(st_bbox(DE_map)), maptype = "stamen_terrain", zoom = 6)
DE_map <-  get_stadiamap(bbox = c(5.5, 47, 15.5, 55), maptype = "stamen_terrain", zoom = 6)

Briwecs <- ggmap(DE_map) +
  geom_point(data = data.frame(lon = coords$Lon, lat = coords$Lat), aes(lon, lat),
             inherit.aes = FALSE, col = color, size = 3) +
  geom_text(data = data.frame(lon = coords$Lon, lat = coords$Lat, loc = coords$Loc),
            aes(label = loc),
            hjust = c(+0.5, +0.5, -0.2, -0.2, -0.2),
            vjust = c(-0.6, -0.6, -0.2, -0.2, +0.8),
            color = color, size = 10) +
  xlab(NULL) + ylab(NULL) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) #+
  #annotate("text", x = 6.1, y = 54.6, label = "A", size = 8)

# ggsave(filename = "plots/locations.png", dpi = 300, width = 20, height = 10, units = "cm")

# DROPS data
coords <- read.csv("DROPS/raw_data/11-Info-Study.tab", sep = "\t")
# coords <- read.csv("DROPS/raw_data/3b-Indices_Env_level.csv")[, c("Experiment", "Lat", "Long")]
coords$Loc <- gsub("(.*)[1-9][1-9]", "\\1", coords$StudyUniqueID)
coords <- unique(coords[, c("Loc", "GeographicLocationLatitude", "GeographicLocationLongitude")])
envs <- rownames(readRDS("DROPS/data/EC.rds"))
envs <- unique(gsub("(.*)[1-9][1-9]", "\\1", envs))
coords <- coords[match(envs, coords$Loc),]
names(coords) <- c("Loc", "Lat", "Lon")
coords.eu <- coords[coords$Loc != "Gra",]

EU_map <- rnaturalearth::ne_countries(continent = "europe", returnclass = "sf")
EU_map <-  get_stadiamap(bbox = c(-2, 40, 30, 52), maptype = "stamen_terrain", zoom = 6)

DROPS <- ggmap(EU_map) +
  geom_point(data = data.frame(lon = coords.eu$Lon, lat = coords.eu$Lat), aes(lon, lat),
             inherit.aes = FALSE, col = color, size = 3) +
  geom_text(data = data.frame(lon = coords.eu$Lon, lat = coords.eu$Lat, loc = coords.eu$Loc),
            aes(label = loc),
            hjust = c(-0.3, +0.5, -0.2, -0.2, +0.1, +0.5, -0.3, -0.3, +0.5),
            vjust = c(+0.5, -0.6, -0.2, -0.2, -0.5, -0.6, +0.5, +0.5, +1.6),
            color = color, size = 10) +
  xlab(NULL) + ylab(NULL) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) #+
  #annotate("text", x = -1, y = 51.42, label = "B", size = 8)

DROPS

cowplot::plot_grid(Briwecs, DROPS, rel_widths = c(1, 2.29)) +
  theme(plot.margin = grid::unit(c(0, 0, 0, 0), "pt"),
        plot.background = element_rect(fill = "white", color = NA))

ggsave(filename = "plots/locations.png", dpi = 300, width = 30, height = 15, units = "cm")

rm(list = ls())

# Plotting the Gaussian kernel function ====
library(ggplot2)
distance <- bigsnpr::seq_log(0.01, 5, 50)
bandwidth <- bigsnpr::seq_log(0.001, 2, 50)
x <- expand.grid(Distance = distance, Bandwidth = bandwidth)
x$Correlation <- exp(-x$Bandwidth * x$Distance)

# Custom color palette:
my_palette <- colorRampPalette(c("#fcdd06", "#db161f", "#0e44af"))
my_colors <- my_palette(10000)

ggplot(x, aes(x = Distance, y = Bandwidth, z = Correlation)) +
  geom_contour_filled(bins = 10000) +
  geom_contour(bins = 10, color = "black", alpha = 0.5) +
  geom_point(data = x, aes(x = Distance, y = Bandwidth, color = Correlation), alpha = 0) +
  theme_classic(base_size = 30) +
  xlab("Distance") + ylab("Bandwidth (h)") +
  labs(fill = "Correlation") +
  scale_fill_manual(values = my_colors[10000:1], guide = "none") +
  scale_color_gradientn(colors = my_colors[10000:1],
                        limits = c(0, 1)) +
  theme(legend.title = element_text(vjust = 5),
        legend.key.width = unit(3, "cm"),
        legend.key.height = unit(1.5, "cm"))

# Export manually via Export > Save as Image at 2000 x 750 (width x height) pixels because the colors are messed up when saving via ggsave()




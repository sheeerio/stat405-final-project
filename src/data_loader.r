required_pkgs <- c("sf", "dplyr", "ggplot2", "readr")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)

library(sf)
library(dplyr)
library(ggplot2)
library(readr)

nfdb <- st_read("data/NFDB_point/NFDB_point_20250519.shp")
coords <- st_coordinates(nfdb)
nfdb_df <- nfdb %>% st_drop_geometry()
nfdb_df$LONGITUDE <- coords[,1]
nfdb_df$LATITUDE <- coords[,2]

nfdb_df <- as.data.frame(nfdb)
nfdb_df <- nfdb_df[!is.na(nfdb_df$SIZE_HA) & nfdb_df$SIZE_HA > 0, ]
# structure
str(nfdb)

# first 10 rows
print(head(nfdb, 10))

# dims
cat("Rows:", nrow(nfdb), " Columns:", ncol(nfdb), "\n")

# col names
print(names(nfdb))

# summary
if ("SIZE_HA" %in% names(nfdb)) {
  cat("\nFire size (hectares) summary:\n")
  print(summary(nfdb$SIZE_HA))
}

if ("YEAR" %in% names(nfdb)) {
  cat("\nYear range:\n")
  cat("Min:", min(nfdb$YEAR, na.rm = TRUE), 
      " Max:", max(nfdb$YEAR, na.rm = TRUE), "\n")
}

if ("CAUSE" %in% names(nfdb)) {
  cat("\nFire causes:\n")
  print(table(nfdb$CAUSE))
}

nfdb_df <- nfdb %>% st_drop_geometry()

if ("SIZE_HA" %in% names(nfdb_df)) {
  p1 <- ggplot(nfdb_df %>% filter(SIZE_HA > 0), aes(x = SIZE_HA)) +
    geom_histogram(bins = 50, fill = "darkorange", color = "black", alpha = 0.7) +
    scale_x_log10() +
    labs(
      title = "Distribution of Wildfire Sizes in Canada (NFDB)",
      x = "Fire Size (hectares, log scale)",
      y = "Count"
    ) +
    theme_minimal()
  
  ggsave("figures/fire_size_distribution.png", p1, width = 8, height = 5, dpi = 150)
  print(p1)
}

if ("YEAR" %in% names(nfdb_df)) {
  fires_per_year <- nfdb_df %>%
    filter(!is.na(YEAR)) %>%
    count(YEAR)
  
  p2 <- ggplot(fires_per_year, aes(x = YEAR, y = n)) +
    geom_line(color = "firebrick", linewidth = 0.8) +
    geom_point(color = "firebrick", size = 1) +
    labs(
      title = "Number of Wildfires per Year in Canada",
      x = "Year",
      y = "Number of Fires"
    ) +
    theme_minimal()
  
  ggsave("figures/fires_per_year.png", p2, width = 8, height = 5, dpi = 150)
  print(p2)
}

if (all(c("LATITUDE", "LONGITUDE") %in% names(nfdb_df))) {
  set.seed(42)
  sample_df <- nfdb_df %>%
    filter(!is.na(LATITUDE), !is.na(LONGITUDE)) %>%
    slice_sample(n = min(5000, nrow(.)))
  
  p3 <- ggplot(sample_df, aes(x = LONGITUDE, y = LATITUDE)) +
    geom_point(alpha = 0.3, size = 0.5, color = "red") +
    labs(
      title = "Spatial Distribution of Wildfires (sample)",
      x = "Longitude",
      y = "Latitude"
    ) +
    coord_quickmap() +
    theme_minimal()
  
  ggsave("figures/fire_locations_map.png", p3, width = 8, height = 5, dpi = 150)
  print(p3)
}
library(sf)
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(scales)
library(patchwork)   # install.packages("patchwork") if needed

dir.create("figures", showWarnings = FALSE)

nfdb <- st_read("data/NFDB_point/NFDB_point_20250519.shp")
coords <- st_coordinates(nfdb)
df <- nfdb %>%
  st_drop_geometry() %>%
  mutate(
    LONGITUDE = coords[, 1],
    LATITUDE  = coords[, 2]
  ) %>%
  filter(!is.na(SIZE_HA), SIZE_HA > 0)

cat("Rows after filtering:", nrow(df), "\n")

# 1a. Density by CAUSE
df <- df %>%
  mutate(CAUSE_LABEL = case_when(
    CAUSE == "H"  ~ "Human",
    CAUSE == "L"  ~ "Lightning",
    CAUSE == "H-PB" ~ "Prescribed Burn",
    TRUE          ~ "Other/Unknown"
  ))

p_cause <- ggplot(df %>% filter(CAUSE_LABEL %in% c("Human", "Lightning")),
                  aes(x = SIZE_HA, fill = CAUSE_LABEL)) +
  geom_density(alpha = 0.5, adjust = 0.8) +
  scale_x_log10(labels = comma) +
  scale_fill_manual(values = c("Human" = "steelblue", "Lightning" = "darkorange")) +
  labs(
    title    = "Fire Size Distribution by Cause",
    subtitle = "Bimodality may partly reflect two populations: human vs. lightning fires",
    x        = "Fire Size (hectares, log scale)",
    y        = "Density",
    fill     = "Cause"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

ggsave("figures/bimodality_by_cause.png", p_cause, width = 8, height = 5, dpi = 150)
print(p_cause)

# 1b. Density by reporting agency (province) — top 8 agencies
top_agencies <- df %>%
  count(SRC_AGENCY, sort = TRUE) %>%
  slice_head(n = 8) %>%
  pull(SRC_AGENCY)

p_agency <- df %>%
  filter(SRC_AGENCY %in% top_agencies) %>%
  ggplot(aes(x = SIZE_HA, colour = SRC_AGENCY)) +
  geom_density(adjust = 0.8, linewidth = 0.8) +
  scale_x_log10(labels = comma) +
  facet_wrap(~SRC_AGENCY, ncol = 4) +
  labs(
    title    = "Fire Size Distribution by Reporting Agency",
    subtitle = "Reporting thresholds vary by province — can produce artificial modes",
    x        = "Fire Size (hectares, log scale)",
    y        = "Density"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

ggsave("figures/bimodality_by_agency.png", p_agency, width = 12, height = 6, dpi = 150)
print(p_agency)

# 1c. Minimum reported size per agency — reveals reporting thresholds
min_size_by_agency <- df %>%
  group_by(SRC_AGENCY) %>%
  summarise(
    min_size  = min(SIZE_HA),
    n_fires   = n(),
    pct_small = mean(SIZE_HA < 1) * 100
  ) %>%
  arrange(min_size)

cat("\nMinimum reported fire size by agency (first 5 rows):\n")
print(head(min_size_by_agency, 5))

annual <- df %>%
  filter(!is.na(YEAR), YEAR >= 1959, YEAR <= 2024) %>%
  group_by(YEAR) %>%
  summarise(
    n_fires        = n(),
    total_area     = sum(SIZE_HA, na.rm = TRUE),
    mean_size      = mean(SIZE_HA, na.rm = TRUE),
    median_size    = median(SIZE_HA, na.rm = TRUE),
    n_large        = sum(SIZE_HA >= 10000, na.rm = TRUE),
    pct_large      = mean(SIZE_HA >= 10000, na.rm = TRUE) * 100,
    mean_lat       = mean(LATITUDE, na.rm = TRUE),
    mean_month     = mean(MONTH, na.rm = TRUE),
    .groups = "drop"
  )

# 2a. Total area burned per year (the headline "worse" metric)
p_area <- ggplot(annual, aes(x = YEAR, y = total_area / 1e6)) +
  geom_col(fill = "firebrick", alpha = 0.7) +
  geom_smooth(method = "lm", colour = "black", linewidth = 0.8, se = TRUE) +
  scale_y_continuous(labels = function(x) paste0(x, "M")) +
  labs(
    title    = "Total Area Burned per Year (Canada)",
    subtitle = "Linear trend overlaid — tests whether aggregate burn area is increasing",
    x        = "Year",
    y        = "Total Area Burned (million hectares)"
  ) +
  theme_minimal(base_size = 13)

ggsave("figures/total_area_per_year.png", p_area, width = 9, height = 5, dpi = 150)
print(p_area)

# 2b. Mean fire size per year (conditional on ignition)
p_mean <- ggplot(annual, aes(x = YEAR, y = mean_size)) +
  geom_line(colour = "darkorange", linewidth = 0.7) +
  geom_smooth(method = "lm", colour = "black", linewidth = 0.8, se = TRUE) +
  scale_y_log10(labels = comma) +
  labs(
    title    = "Mean Fire Size per Year (log scale)",
    subtitle = "Is each fire getting bigger on average?",
    x        = "Year",
    y        = "Mean Fire Size (ha, log scale)"
  ) +
  theme_minimal(base_size = 13)

ggsave("figures/mean_size_per_year.png", p_mean, width = 9, height = 5, dpi = 150)
print(p_mean)

# 2c. Number of large fires (>10,000 ha) per year
p_large <- ggplot(annual, aes(x = YEAR, y = n_large)) +
  geom_col(fill = "darkred", alpha = 0.8) +
  geom_smooth(method = "lm", colour = "black", linewidth = 0.8, se = TRUE) +
  labs(
    title    = "Number of Large Fires (>10,000 ha) per Year",
    subtitle = "Extreme fire frequency as a definition of 'worse'",
    x        = "Year",
    y        = "Count of Large Fires"
  ) +
  theme_minimal(base_size = 13)

ggsave("figures/large_fires_per_year.png", p_large, width = 9, height = 5, dpi = 150)
print(p_large)

lm_area   <- lm(total_area   ~ YEAR, data = annual)
lm_mean   <- lm(log(mean_size) ~ YEAR, data = annual)
lm_large  <- lm(n_large      ~ YEAR, data = annual)

cat("Total area ~ Year: slope =", round(coef(lm_area)[2], 0),
    "ha/yr, p =", round(summary(lm_area)$coefficients[2, 4], 4), "\n")
cat("log(mean size) ~ Year: slope =", round(coef(lm_mean)[2], 5),
    ", p =", round(summary(lm_mean)$coefficients[2, 4], 4), "\n")
cat("Large fires ~ Year: slope =", round(coef(lm_large)[2], 3),
    ", p =", round(summary(lm_large)$coefficients[2, 4], 4), "\n")

# 3a. Mean latitude per year (northward shift?)
p_lat <- ggplot(annual, aes(x = YEAR, y = mean_lat)) +
  geom_line(colour = "steelblue", linewidth = 0.7) +
  geom_smooth(method = "lm", colour = "black", linewidth = 0.8, se = TRUE) +
  labs(
    title    = "Mean Latitude of Wildfires per Year",
    subtitle = "Are fires shifting poleward over time?",
    x        = "Year",
    y        = "Mean Latitude (°N)"
  ) +
  theme_minimal(base_size = 13)

ggsave("figures/mean_latitude_per_year.png", p_lat, width = 9, height = 5, dpi = 150)
print(p_lat)

lm_lat <- lm(mean_lat ~ YEAR, data = annual)
cat("\nLatitude ~ Year: slope =", round(coef(lm_lat)[2], 5),
    "°/yr,  p =", round(summary(lm_lat)$coefficients[2, 4], 4), "\n")

# 3b. Seasonal timing — mean month of fire per year (earlier season?)
p_month <- ggplot(annual, aes(x = YEAR, y = mean_month)) +
  geom_line(colour = "forestgreen", linewidth = 0.7) +
  geom_smooth(method = "lm", colour = "black", linewidth = 0.8, se = TRUE) +
  scale_y_continuous(breaks = 1:12,
                     labels = month.abb) +
  labs(
    title    = "Mean Month of Wildfire per Year",
    subtitle = "Is the fire season starting earlier?",
    x        = "Year",
    y        = "Mean Fire Month"
  ) +
  theme_minimal(base_size = 13)

ggsave("figures/mean_month_per_year.png", p_month, width = 9, height = 5, dpi = 150)
print(p_month)

cause_annual <- df %>%
  filter(!is.na(YEAR), YEAR >= 1970, YEAR <= 2024,
         CAUSE_LABEL %in% c("Human", "Lightning")) %>%
  count(YEAR, CAUSE_LABEL) %>%
  group_by(YEAR) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

p_cause_trend <- ggplot(cause_annual, aes(x = YEAR, y = pct, colour = CAUSE_LABEL)) +
  geom_line(linewidth = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, linetype = "dashed") +
  scale_colour_manual(values = c("Human" = "steelblue", "Lightning" = "darkorange")) +
  labs(
    title    = "Cause Composition of Fires Over Time",
    subtitle = "% of fires attributed to human vs. lightning ignition per year",
    x        = "Year",
    y        = "Percentage of Fires (%)",
    colour   = "Cause"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

ggsave("figures/cause_trend_over_time.png", p_cause_trend, width = 9, height = 5, dpi = 150)
print(p_cause_trend)

# 3d. Decade-level area burned by cause (lightning vs. human — who burns more?)
decade_cause <- df %>%
  filter(!is.na(YEAR), YEAR >= 1970,
         CAUSE_LABEL %in% c("Human", "Lightning")) %>%
  mutate(DECADE = paste0(floor(YEAR / 10) * 10, "s")) %>%
  group_by(DECADE, CAUSE_LABEL) %>%
  summarise(total_area = sum(SIZE_HA) / 1e6, .groups = "drop")

p_decade <- ggplot(decade_cause, aes(x = DECADE, y = total_area, fill = CAUSE_LABEL)) +
  geom_col(position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = c("Human" = "steelblue", "Lightning" = "darkorange")) +
  labs(
    title    = "Total Area Burned by Decade and Cause",
    x        = "Decade",
    y        = "Total Area (million ha)",
    fill     = "Cause"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

ggsave("figures/decade_area_by_cause.png", p_decade, width = 9, height = 5, dpi = 150)
print(p_decade)

summary_by_decade <- df %>%
  filter(!is.na(YEAR), YEAR >= 1970, YEAR <= 2024) %>%
  mutate(DECADE = paste0(floor(YEAR / 10) * 10, "s")) %>%
  group_by(DECADE) %>%
  summarise(
    n_fires       = n(),
    total_area_Mha = round(sum(SIZE_HA) / 1e6, 2),
    mean_size_ha   = round(mean(SIZE_HA), 0),
    n_large        = sum(SIZE_HA >= 10000),
    mean_lat       = round(mean(LATITUDE, na.rm = TRUE), 2),
    .groups = "drop"
  )

print(summary_by_decade)

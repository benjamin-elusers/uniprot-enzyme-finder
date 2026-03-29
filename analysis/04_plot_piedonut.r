library(dplyr)
library(ggplot2)
library(ggrepel)

# 1. Helper function
format_k_m <- function(x) {
  case_when(x >= 1e6 ~ sprintf("%.1f M", x / 1e6), x >= 1e3 ~ sprintf("%.1f K", x / 1e3), TRUE ~ as.character(x))
}

# 2. Load data
dat <- readRDS("data/plot_data/aggregated_plot_data.rds")

# 3. Hierarchical Palette
combined_palette <- c(
  # --- INNER RING BASES ---
  "Full EC"    = "#3182bd",
  "Partial EC" = "#31a354",
  "No EC"      = "#de2d26",

  # --- OUTER RING: FULL EC (Blues) ---
  "Full EC - Enzyme"           = "#3182bd",
  "Full EC - Suspected_Enzyme" = "#6baed6",
  "Full EC - Pseudo_Enzyme"    = "#9ecae1",
  "Full EC - Non_Enzyme"       = "#c6dbef",

  # --- OUTER RING: PARTIAL EC (Greens) ---
  "Partial EC - Enzyme"           = "#31a354",
  "Partial EC - Suspected_Enzyme" = "#74c476",
  "Partial EC - Pseudo_Enzyme"    = "#a1d99b",
  "Partial EC - Non_Enzyme"       = "#c7e9c0",

  # --- OUTER RING: NO EC (Reds) ---
  "No EC - Enzyme"           = "#de2d26",
  "No EC - Suspected_Enzyme" = "#fb6a4a",
  "No EC - Pseudo_Enzyme"    = "#fc9272",
  "No EC - Non_Enzyme"       = "#fcbba1"
)

# 4. MATHEMATICAL LAYOUT PREPARATION

# --- Define Radii, Explosion, and Padding ---
hole_radius <- 3.0
inner_width <- 1.5
gap_width   <- 0.25
outer_width <- 2.0

explode_target <- c("Full EC","Partial EC")
explode_shift  <- 0.3    # Radial explosion offset
angular_pad    <- 0.0035 # Angular gap to detach slices (0.15% of the circle)

# A. Outer Ring (Primary Class)
df_outer <- dat$master %>%
  group_by(ec_cat, primary_class) %>%
  summarize(count = sum(count), .groups = "drop") %>%
  mutate(
    ec_cat = factor(ec_cat, levels = c("Full EC", "Partial EC", "No EC")),
    primary_class = factor(primary_class, levels = c("Enzyme", "Suspected_Enzyme", "Pseudo_Enzyme", "Non_Enzyme"))
  ) %>%
  arrange(ec_cat, primary_class) %>%
  mutate(
    fraction = count / sum(count),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    mid_y = (ymin + ymax) / 2,

    # -----------------------------------------------------------------
    # pmin() ensures we never invert a slice that is smaller than the gap itself
    plot_ymin = ymin + pmin(angular_pad, fraction * 0.4),
    plot_ymax = ymax - pmin(angular_pad, fraction * 0.4),
    stroke_w  = if_else(fraction < 0.005, 0, 0.5), # Thin borders for <0.5% chunks
    # -----------------------------------------------------------------

    explode = if_else(ec_cat  %in%  explode_target, explode_shift, 0),
    xmin = hole_radius + inner_width + gap_width + explode,
    xmax = hole_radius + inner_width + gap_width + outer_width + explode,
    mid_x = (xmin + xmax) / 2,

    pct = fraction * 100,
    pct_text = if_else(pct > 1.0, sprintf("%.1f%%", pct), ""),
    label_text = if_else(pct > 1.0, as.character(primary_class), ""),
    fill_category = paste(ec_cat, primary_class, sep = " - ")
  )

# B. Inner Ring (EC Category)
df_inner <- df_outer %>%
  group_by(ec_cat) %>%
  summarize(
    ymin = min(ymin),
    ymax = max(ymax),
    count = sum(count),
    explode = first(explode),
    .groups = "drop"
  ) %>%
  mutate(
    mid_y = (ymin + ymax) / 2,
    fraction = count / sum(count),

    # Apply the exact same angular padding to the inner ring
    plot_ymin = ymin + pmin(angular_pad, fraction * 0.4),
    plot_ymax = ymax - pmin(angular_pad, fraction * 0.4),

    xmin = hole_radius + explode,
    xmax = hole_radius + inner_width + explode,
    mid_x = (xmin + xmax) / 2,

    pct = fraction * 100,
    pct_text = if_else(pct > 2.0, sprintf("%.1f%%", pct), ""),
    label_text = as.character(ec_cat),
    fill_category = as.character(ec_cat)
  )

# 5. PLOT GENERATION
g_native_donut <- ggplot() +

  # --- GEOMETRIES ---
  # Outer Ring (Using plot_ymin/ymax for gaps, and dynamic linewidth)
  geom_rect(data = df_outer,
            aes(ymin = plot_ymin, ymax = plot_ymax, xmin = xmin, xmax = xmax, fill = fill_category, linewidth = stroke_w),
            color = "white") +

  # Inner Ring
  geom_rect(data = df_inner,
            aes(ymin = plot_ymin, ymax = plot_ymax, xmin = xmin, xmax = xmax, fill = fill_category),
            color = "white", linewidth = 0.8) +

  scale_linewidth_identity() + # Tells ggplot to use the exact stroke_w numbers we calculated

  # --- TEXT AND LABELS ---
  # 0. Center Label (N = 203.1M entries)
  annotate("text", x = 0, y = 0,
           label = sprintf("N = %s\nentries", format_k_m(dat$n_total)),
           size = 18 / .pt, fontface = "bold", lineheight = 0.9,
           hjust = 0.5, vjust = 0.5, color = "black") +

  # 1. Inner Percentages
  geom_text(data = df_inner %>% filter(pct_text != ""),
            aes(x = mid_x, y = mid_y, label = pct_text),
            size = 14 / .pt, fontface = "bold", color = "white") +

  # 2. Outer Percentages
  geom_text(data = df_outer %>% filter(pct_text != ""),
            aes(x = mid_x, y = mid_y, label = pct_text),
            size = 12 / .pt, fontface = "bold", color = "grey20") +

  # 3. Inner Labels
  geom_text(data = df_inner,
            aes(x = xmin - 1.2, y = mid_y, label = label_text, color = fill_category),
            size = 18 / .pt, fontface = "italic", show.legend = FALSE) +

  # 4. Outer Labels
  geom_text_repel(data = df_outer %>% filter(label_text != ""),
                  aes(x = xmax + 0.2, y = mid_y, label = label_text),
                  color = "grey20", size = 14 / .pt, fontface = "bold",
                  segment.color = "grey50", direction = "y", nudge_x = 0.5, show.legend = FALSE) +

  # --- COORDINATES & THEMING ---
  coord_polar(theta = "y", direction = -1) +
  xlim(c(0, 8)) +

  scale_fill_manual(values = combined_palette) +
  scale_color_manual(values = combined_palette) +

  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 16, hjust = 0.5, margin = margin(b = 15))
  ) +
  labs(
    title = "UniProt Enzyme Landscape",
    subtitle = "Uniprot release 2026_01"
  )

print(g_native_donut)

# Export
dir.create("plots", showWarnings = FALSE)
ggsave("plots/04_enzyme_ec_donut.pdf", g_native_donut, width = 8, height = 8, bg = "white")
ggsave("plots/04_enzyme_ec_donut.png", g_native_donut, width = 8, height = 8, bg = "white")


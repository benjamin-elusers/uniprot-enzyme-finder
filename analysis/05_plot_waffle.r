library(dplyr)
library(tidyr)
library(ggplot2)

format_k_m <- function(x) { case_when(x >= 1e6 ~ sprintf("%.1f M", x / 1e6), x >= 1e3 ~ sprintf("%.1f K", x / 1e3), TRUE ~ as.character(x)) }

dat <- readRDS("data/plot_data/aggregated_plot_data.rds")
blue_palette <- c("Full EC" = "#08519c", "Partial EC" = "#3182bd", "No EC" = "#9ecae1")
class_palette <- c("Enzyme" = "#44cf32", "Suspected_Enzyme" = "#a185f8", "Pseudo_Enzyme" = "#df2f3f", "Non_Enzyme" = "#c9c9c9")

sq_val <- 500000
df_waffle <- dat$master %>%
  mutate(
    has_ec = factor(if_else(ec_cat == "No EC", "No EC", "With EC"), levels=c("With EC", "No EC")),
    class = case_when(
      primary_class %in% c("Enzyme", "Pseudo_Enzyme") ~ "Enzyme",
      primary_class == "Suspected_Enzyme" ~ "Suspected_Enzyme",
      TRUE ~ "Non_Enzyme"
    ),
    squares = ceiling(round(count / sq_val))
  ) %>%
  filter(squares > 0) %>% uncount(squares) %>%
  group_by(has_ec) %>%
  mutate(square_id = row_number() - 1, x = square_id %% 10, y = square_id %/% 10) %>%
  ungroup() %>%
  mutate(fill_category = if_else(has_ec == "With EC", if_else(ec_cat == "Full EC", "Full EC", "Partial EC"), class))

g_waffle <- ggplot(df_waffle, aes(x = x, y = y, fill = primary_class, color = ec_cat)) +
  geom_tile(linewidth = 1, alpha = 0.8) +
  facet_wrap(~ has_ec, scales = "free_x") +
  scale_fill_manual(values = class_palette) +
  scale_color_manual(values = blue_palette, na.value = "white") +
  guides(
    fill = guide_legend(nrow = 2, byrow = TRUE),
    color = guide_legend(nrow = 1, override.aes = list(fill = "white", linewidth = 1))
  ) +
  theme_void() +
  theme(
    strip.text = element_text(size = 14, face = "bold", margin = margin(b = 10)),
    legend.position = "top", legend.box = "vertical", legend.margin = margin(b = 10),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15))
  ) +
  labs(title = "Enzyme Landscape", fill = "Category", color = "EC Status")

print(g_waffle)
ggsave("plots/05_enzyme_ec_waffle.png", g_waffle, width=6, height=6, bg='white')
ggsave("plots/05_enzyme_ec_waffle.pdf", g_waffle, width=6, height=6, bg='white')


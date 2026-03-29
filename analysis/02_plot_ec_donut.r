library(dplyr)
library(ggplot2)
library(ggrepel)

format_k_m <- function(x) { case_when(x >= 1e6 ~ sprintf("%.1f M", x / 1e6), x >= 1e3 ~ sprintf("%.1f K", x / 1e3), TRUE ~ as.character(x)) }

dat <- readRDS("data/plot_data/aggregated_plot_data.rds")
blue_palette <- c("Full EC" = "#08519c", "Partial EC" = "#3182bd", "No EC" = "#9ecae1")

df_ec <- dat$master %>%
  group_by(ec_cat) %>%
  summarize(count = sum(count), .groups="drop") %>%
  mutate(
    ec_cat = factor(ec_cat, levels = c("No EC", "Partial EC", "Full EC")),
    pct = count / sum(count) * 100,
    label_text = paste0(round(pct, 1), "%\n", format_k_m(count))
  ) %>%
  arrange(desc(ec_cat)) %>%
  mutate(csum = cumsum(count), mid_y = csum - (count / 2))

g_ec_donut <- ggplot(df_ec, aes(x = 2.5, y = count, fill = ec_cat)) +
  geom_col(color = "white", linewidth = 0.5) +
  geom_text_repel(aes(x = 3.5, y = mid_y, label = label_text, color = ec_cat),
                  show.legend = FALSE, size = 4, fontface = "bold", segment.color = "grey70", direction = "y", force = 2) +
  coord_polar(theta = "y", direction = -1) +
  xlim(c(0, 3.5)) +
  scale_fill_manual(values = blue_palette, aesthetics = c('fill','color')) +
  theme_void() +
  theme(legend.position = "bottom", legend.title.position = 'top', legend.title = element_text(hjust=0.5, face='bold')) +
  labs(fill = "EC?") +
  annotate("text", x = 0, y = 0, label = paste0("N = ", format_k_m(sum(df_ec$count))), size = 5, fontface = "bold")

print(g_ec_donut)
ggsave("plots/02_ec_donut.png", g_ec_donut, width=6, height=6, bg='white')
ggsave("plots/02_ec_donut.pdf", g_ec_donut, width=6, height=6, bg='white')


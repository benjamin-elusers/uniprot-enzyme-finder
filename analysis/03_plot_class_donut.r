library(dplyr)
library(ggplot2)
library(ggrepel)

format_k_m <- function(x) { case_when(x >= 1e6 ~ sprintf("%.1f M", x / 1e6), x >= 1e3 ~ sprintf("%.1f K", x / 1e3), TRUE ~ as.character(x)) }

dat <- readRDS("data/plot_data/aggregated_plot_data.rds")
class_palette <- c("Enzyme" = "#44cf32", "Suspected_Enzyme" = "#a185f8", "Pseudo_Enzyme" = "#df2f3f", "Non_Enzyme" = "#c9c9c9")

df_class <- dat$master %>%
  group_by(primary_class) %>%
  summarize(count = sum(count), .groups = "drop") %>%
  arrange(count) %>%
  mutate(
    primary_class = factor(primary_class, levels = primary_class),
    pct = count / sum(count) * 100,
    label_text = paste0(round(pct, 1), "%\n", format_k_m(count))
  ) %>%
  arrange(desc(primary_class)) %>%
  mutate(csum = cumsum(count), mid_y = csum - (count / 2))

g_class_donut <- ggplot(df_class, aes(x = 2.5, y = count, fill = primary_class)) +
  geom_col(color = "white", linewidth = 0.2) +
  geom_text_repel(aes(x = 3.5, y = mid_y, label = label_text, color = primary_class),
                  show.legend = FALSE, size = 4, fontface = "bold", segment.color = "grey70", direction = "y", force = 2) +
  coord_polar(theta = "y", direction = 1) +
  xlim(c(0, 3.5)) +
  scale_fill_manual(values = class_palette, aesthetics = c('color','fill')) +
  theme_void() +
  theme(legend.position = "bottom", legend.title.position = 'top', legend.title = element_text(hjust=0.5, face='bold')) +
  labs(fill = "Classification") +
  annotate("text", x = 0, y = 0, label = paste0("N = ", format_k_m(sum(df_class$count))), size = 5, fontface = "bold")

print(g_class_donut)
ggsave("plots/03_enzyme_donut.png", g_class_donut, width=6, height=6, bg='white')
ggsave("plots/03_enzyme_donut.pdf", g_class_donut, width=6, height=6, bg='white')

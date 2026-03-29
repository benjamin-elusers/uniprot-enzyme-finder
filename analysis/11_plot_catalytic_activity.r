library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)     # Required for separate_rows()
library(tidytext)  # Required for reorder_within()

# 1. Helper function
format_k_m <- function(x) {
  case_when(
    x >= 1e6 ~ sprintf("%.1f M", x / 1e6),
    x >= 1e3 ~ sprintf("%.1f K", x / 1e3),
    TRUE ~ as.character(x)
  )
}


# Themes
# Sizes & Layout
THICKNESS_MIN   <- 2.5   # Inner Radius (Hole) -> Decreased to make donut THICKER
THICKNESS_MAX   <- 4.5   # Outer Radius (Ring Edge)
LABEL_OUT_POS   <- 5.5   # Position for outside labels
LIMIT_MAX       <- 6.0   # Plot limit to prevent clipping
SMART_THRESHOLD <- 0.04  # 4% threshold to move labels outside

# Text Sizes (Harmonized)
TOTAL_TEXT_SIZE <- 5
DONUT_TEXT_SIZE <- 3.5
BAR_TEXT_SIZE   <- 3.5
TITLE_SIZE      <- 11
SUBTITLE_SIZE   <- 9
my_theme <- theme_minimal(base_size = 10) +
  theme(
    plot.title    = element_text(face = "bold", size = TITLE_SIZE, color = "#2c3e50", hjust = 0),
    plot.subtitle = element_text(size = SUBTITLE_SIZE, color = "#7f8c8d"),
    legend.position = "none",
    panel.grid.major = element_line(color = "grey85", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.line  = element_line(color = "black", linewidth = 0.2),
    axis.text  = element_text(color = "black", size = 8),
    axis.title = element_blank(),
    plot.margin = margin(3, 3, 3, 3, "mm") # Tight margins
  )
# 2. Load data
dat <- readRDS("data/plot_data/aggregated_plot_data.rds")


p_cat = dat$cat_act %>%
  mutate(has_catalytic_activity = factor(has_catalytic_activity, levels=c(F,T),labels=c('No','Yes'))) %>%
  filter(is_enzyme) %>%
  ungroup() %>%
  dplyr::select(has_catalytic_activity,has_ec,count_by_ec) %>% distinct() %>%
  group_by(has_ec) %>%
  mutate(percentage = count_by_ec / sum(count_by_ec),
         label_count = format_k_m(count_by_ec),
         label_perc = scales::percent(percentage, 0.1)
  ) %>%
  ggplot(aes(fill=has_catalytic_activity, x=count_by_ec, y=has_ec)) +
  scale_x_continuous(labels=format_k_m, expand=expansion(mult=c(0, 0.2))) +
  geom_col(width=0.5) +
  scale_fill_manual(values=c('Yes'="gray30", 'No'="gray90"))+
  #scale_color_manual(values=c("TRUE"="white", "FALSE"="black"))+
  geom_text(aes(label=label_count, group=has_catalytic_activity), hjust=-0.1, size=5, position = position_stack(vjust = 1)) +
  geom_text(aes(label=label_perc, x=cumsum(count_by_ec)), vjust=-1.5, hjust=-0.1, size=5,position = position_stack(vjust =  0)) +
  my_theme +
  theme(panel.grid.major.y=element_blank(),
        axis.line.y = element_blank(),
        aspect.ratio = 0.25,
        axis.text = element_text(size=unit(12,'pt')),
        legend.position = 'top',
        legend.title.position = 'top',
        legend.text = element_text(size=unit(10,'pt')),
        legend.title=element_text(hjust=1, size=unit(12,'pt'))) +
  coord_cartesian(clip='off') +
  theme() +
  guides( fill = guide_legend(title='Catalytic Activity'))
p_cat
ggsave("plots/11_catalytic_activity.png", p_cat, width=8, height=6, bg='white')
ggsave("plots/11_catalytic_activity.pdf", p_cat, width=8, height=6, bg='white')

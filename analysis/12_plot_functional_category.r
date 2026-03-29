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

# 2. Load data
dat <- readRDS("data/plot_data/aggregated_plot_data.rds")

# === P3: Functional Category (Bar) ===
func_palette <- c("Unclassified_Enzyme"="#DDC191", "Core_Metabolism"="#DB6D54",
                  "Info_Processing"="#5D9AD1", "Secondary_Metabolism"="#8FC567",
                  "Signaling"="#F1CD40", "Inactive_Non_Enzyme"="#efa99f", "Transport"="#9D82BF")
# Data Prep for Bars
prepare_bar_data <- function(df, count = 'count') {
  library(scales)
  df %>%
    mutate(
      percentage = !!count / sum(!!count),
      label_text = paste0(format_k_m(!!count), " (", percent(percentage, 0.1), ")")
    )
}

#df_func <- enframe(dat$func_cat, name="category", value="count") %>%
#  unnest(count) %>% prepare_bar_data()


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

# Themes
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

df_func = dat$func_cat %>%
  filter(is_enzyme) %>%
  ungroup() %>%
  dplyr::select(functional_category,count_both) %>% distinct() %>%
  group_by(functional_category) %>%
  summarize(count = sum(count_both)) %>%
  prepare_bar_data(count = .[['count']])

p3 <- ggplot(df_func, aes(x=reorder(functional_category, count), y=count, fill=functional_category)) +
  geom_col(width=0.3) + coord_flip(clip='off') +
  scale_fill_manual(values=func_palette) +
  scale_y_continuous(labels=format_k_m, expand=expansion(mult=c(0, 0.4))) +
  geom_text(aes(label=label_text), hjust=-0.1, size=BAR_TEXT_SIZE) +
  my_theme + theme(panel.grid.major.y=element_blank(), axis.text =  element_text(size=unit(12,'pt')), aspect.ratio = 0.7) +
  labs(title="Functional Classification")
print(p3)

ggsave("plots/12_functional_category.png", p3, width=8, height=6, bg='white')
ggsave("plots/12_functional_category.pdf", p3, width=8, height=6, bg='white')


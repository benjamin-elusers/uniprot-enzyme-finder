library(dplyr)
library(stringr)
library(ggplot2)

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

# 3. Abbreviation dictionary
abbreviations <- c(
  "GO_Catalytic"="GO", "Keyword_Catalytic"="KW", "Protein_Name"="Name",
  "Active_Site"="Site", "Reaction_Text"="Rxn", "Rhea_ID"="Rhea",
  "Enzyme_DB"="DB", "Text_Mining"="Text", "EC_Number"="EC", "\\|"=" + "
)

# 4. Data processing
df_plot <- dat$evidence %>%
  mutate(
    ev_count = if_else(evidence_field %in% c("None", ""), 0L, str_count(evidence_field, "\\|") + 1L),
    ev_group = case_when(
      ev_count == 0 ~ "0 Ev",
      ev_count == 1 ~ "1 Evidence",
      ev_count == 2 ~ "2 Evidences",
      TRUE ~ ">2 Evidences"
    ),
    evidence_clean = str_replace_all(evidence_field, abbreviations)
  ) %>%
  filter(ev_group != "0 Ev") %>%
  group_by(has_ec, ev_group) %>%
  mutate(
    total_in_bucket = sum(count),
    pct_of_bucket = (count / total_in_bucket) * 100,
    plot_category = if_else(pct_of_bucket >= 2.0, evidence_clean, "Other (<2%)")
  ) %>%
  group_by(has_ec, ev_group, plot_category, total_in_bucket) %>%
  summarize(count = sum(count), .groups = "drop") %>%
  mutate(
    pct_final = (count / total_in_bucket) * 100,
    # Combine the abbreviation and the percentage for the text label
    label_text = if_else(pct_final > 4.0, paste0(plot_category, "\n", round(pct_final, 1), "%"), ""),
    has_ec = factor(has_ec, levels = c("With EC", "No EC")),
    ev_group = factor(ev_group, levels = c("1 Evidence", "2 Evidences", ">2 Evidences")),
    plot_category = str_remove_all(string = plot_category, pattern = "EC \\+ "),
    # Apply Confidence Tiers
    confidence_tier = case_when(
      plot_category == "EC" ~ "1. Enzyme activity",
      str_detect(plot_category, "Site|Rxn|Rhea|DB") ~ "2. High (Biochemical/Structural)",
      str_detect(plot_category, "GO|KW") ~ "3. Medium (Ontology/Homology)",
      str_detect(plot_category, "Name|Text") ~ "4. Low (Lexical/Text)",
      TRUE ~ "5. Other/Mixed"
    )
  ) %>%
  # Arrange so that chunks of the same confidence tier stack neatly next to each other
  arrange(has_ec, ev_group, desc(confidence_tier), count) %>%
  mutate(plot_category = factor(plot_category, levels = unique(plot_category)))

# Create dataset for the absolute totals at the end of the bars
df_totals <- df_plot %>% group_by(has_ec, ev_group, total_in_bucket) %>% slice(1)

# 5. Semantic Palette
tier_palette <- c(
  "1. Enzyme activity"               = "#cc1616", # Red Enzyme
  "2. High (Biochemical/Structural)" = "#2b83ba", # Strong Blue
  "3. Medium (Ontology/Homology)"    = "#abdda4", # Soft Green
  "4. Low (Lexical/Text)"            = "#fdae61", # Warning Orange
  "5. Other/Mixed"                   = "#cccccc"  # Neutral Grey
)

# 6. Plotting
g_stacked <- ggplot(df_plot, aes(x = count, y = has_ec, fill = confidence_tier, group = plot_category)) +
  # Use white lines to clearly separate individual combinations inside the same color block
  geom_col(position = "stack", color = "white", linewidth = 0.3) +

  # Text inside the bars (changed to grey10 so it's readable against the light green/orange backgrounds)
  geom_text(aes(label = label_text), position = position_stack(vjust = 0.5), size = 3, fontface = "bold", color = "grey10") +

  # Total absolute count label at the end of each bar
  geom_text(data = df_totals, aes(x = total_in_bucket, y = has_ec, label = format_k_m(total_in_bucket)),
            inherit.aes = FALSE, hjust = -0.1, size = 3.5, fontface = "bold") +

  # Faceting
  facet_wrap(~ ev_group, ncol = 3, scales = "free_x") +

  # Ensure the total labels don't get cut off
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  scale_fill_manual(values = tier_palette) +

  theme_minimal() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face="bold")
  ) +
  labs(
    title = "Evidence Breakdown by Number of Annotations",
    subtitle = "Bar width is proportional to total entries. Color represents the strongest evidence class present.",
    fill = "Highest Evidence Class"
  )

print(g_stacked)
ggsave("plots/06_evidence_stacked.png", g_stacked, width=14, height=10, bg='white')
ggsave("plots/06_evidence_stacked.pdf", g_stacked, width=14, height=10, bg='white')


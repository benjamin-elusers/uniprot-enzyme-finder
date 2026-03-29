# 08_plot_evidence_tree.R

library(dplyr)
library(stringr)
library(ggplot2)
library(treemapify)


# 1. Helper function for formatting numbers
format_k_m <- function(x) {
  case_when(
    x >= 1e6 ~ sprintf("%.1f M", x / 1e6),
    x >= 1e3 ~ sprintf("%.1f K", x / 1e3),
    TRUE ~ as.character(x)
  )
}

# 2. Load the aggregated master data
dat <- readRDS("data/plot_data/aggregated_plot_data.rds")

# 3. Define the abbreviation dictionary
# Maps the raw python evidence strings to short labels for the treemap
abbreviations <- c(
  "GO_Catalytic"      = "GO",
  "Keyword_Catalytic" = "KW",
  "Protein_Name"      = "Name",
  "Active_Site"       = "Site",
  "Reaction_Text"     = "Rxn",
  "Rhea_ID"           = "Rhea",
  "Enzyme_DB"         = "DB",
  "Text_Mining"       = "Text",
  "EC_Number"         = "EC",
  "\\|"               = " + "  # Replaces the python pipe separator with a plus sign
)


abbreviations <- c("GO_Catalytic"="GO", "Keyword_Catalytic"="KW", "Protein_Name"="Name",
                   "Active_Site"="Site", "Reaction_Text"="Rxn", "Rhea_ID"="Rhea",
                   "Enzyme_DB"="DB", "Text_Mining"="Text", "EC_Number"="EC", "\\|"=" + ")

df_plot <- dat$evidence %>%
  mutate(
    ev_count = if_else(evidence_field %in% c("None", ""), 0L, str_count(evidence_field, "\\|") + 1L),
    ev_group = case_when(ev_count == 0 ~ "0 Ev", ev_count == 1 ~ "1 Evidence",
                         ev_count == 2 ~ "2 Evidences", TRUE ~ ">2 Evidences"),
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
    has_ec = factor(has_ec, levels = c("With EC", "No EC")),
    ev_group = factor(ev_group, levels = c("1 Evidence", "2 Evidences", ">2 Evidences")),
    plot_category = str_remove_all(string = plot_category, pattern = "EC \\+ "),

    confidence_tier = case_when(
      plot_category == "EC" ~ "1. Enzyme activity",
      str_detect(plot_category, "Site|Rxn|Rhea|DB|EC") ~ "2. High (Biochemical/Structural)",
      str_detect(plot_category, "GO|KW") ~ "3. Medium (Ontology/Homology)",
      str_detect(plot_category, "Name|Text") ~ "4. Low (Lexical/Text)",
      TRUE ~ "4. Other/Mixed"
    ),
    tree_label = if_else(pct_final > 3.0, paste0(plot_category, "\n", round(pct_final, 1), "%\n(", format_k_m(count), ")"), "")
  )


# 5. Semantic Palette
tier_palette <- c(
  "1. Enzyme activity"               = "#cc1616", # Red Enzyme
  "2. High (Biochemical/Structural)" = "#2b83ba", # Strong Blue
  "3. Medium (Ontology/Homology)"    = "#abdda4", # Soft Green
  "4. Low (Lexical/Text)"            = "#fdae61", # Warning Orange
  "5. Other/Mixed"                   = "#cccccc"  # Neutral Grey
)
# Plot the Treemap
g_treemap <- ggplot(df_plot, aes(
  area = count,
  fill = confidence_tier,  # <--- Map fill to the new semantic tier!
  subgroup = ev_group,
  label = tree_label
)) +
  geom_treemap(color = "white", size = 0.5) +
  geom_treemap_subgroup_border(color = "black", size = 2) +
  geom_treemap_subgroup_text(place = "bottom", grow = FALSE, alpha = 0.4, colour = "black", fontface = "italic", size = 14) +
  geom_treemap_text(colour = "grey10", place = "centre", reflow = TRUE, fontface = "bold", min.size = 8) +
  facet_wrap(~ has_ec) +

  # Use the custom semantic palette
  scale_fill_manual(values = tier_palette) +

  theme_minimal() +
  theme(
    strip.text = element_text(size = 16, face = "bold", margin = margin(b = 10)),
    # We turn the legend back ON so readers know what the colors mean
    legend.position = "bottom",
    legend.title = element_text(face="bold"),
    plot.title = element_text(size = 18, face = "bold", margin = margin(b = 15)),
    plot.subtitle = element_text(size = 12, margin = margin(b = 15))
  ) +
  labs(
    title = "Evidence Breakdown by Number of Annotations",
    subtitle = "Area is proportional to entry count. Color represents the strongest evidence class present.",
    fill = "Highest Evidence Class"
  )

print(g_treemap)

# ==============================================================================
# TREE 1: Enzymes and Suspected Enzymes WITH EC
# ==============================================================================
df_with_ec <- df_plot %>% filter(has_ec == "With EC")

g_treemap_with_ec <- ggplot(df_with_ec, aes(
  area = count,
  fill = confidence_tier,
  subgroup = ev_group,
  label = tree_label
)) +
  geom_treemap(color = "white", size = 0.5) +
  geom_treemap_subgroup_border(color = "black", size = 2) +
  geom_treemap_subgroup_text(place = "bottom", grow = FALSE, alpha = 0.4, colour = "black", fontface = "italic", size = 14) +
  geom_treemap_text(colour = "grey10", place = "centre", reflow = TRUE, fontface = "bold", min.size = 8) +

  scale_fill_manual(values = tier_palette) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face="bold"),
    plot.title = element_text(size = 18, face = "bold", margin = margin(b = 10)),
    plot.subtitle = element_text(size = 12, margin = margin(b = 15))
  ) +
  labs(
    title = "Evidence Breakdown: Enzymes WITH an EC Number",
    subtitle = "Area is proportional to entry count. Color represents the strongest evidence class present.",
    fill = "Highest Evidence Class"
  )

print(g_treemap_with_ec)

# ==============================================================================
# TREE 2: Enzymes and Suspected Enzymes WITHOUT EC (The "Dark Matter")
# ==============================================================================
df_no_ec <- df_plot %>% filter(has_ec == "No EC")

g_treemap_no_ec <- ggplot(df_no_ec, aes(
  area = count,
  fill = confidence_tier,
  subgroup = ev_group,
  label = tree_label
)) +
  geom_treemap(color = "white", size = 0.5) +
  geom_treemap_subgroup_border(color = "black", size = 2) +
  geom_treemap_subgroup_text(place = "bottom", grow = FALSE, alpha = 0.4, colour = "black", fontface = "italic", size = 14) +
  geom_treemap_text(colour = "grey10", place = "centre", reflow = TRUE, fontface = "bold", min.size = 8) +

  scale_fill_manual(values = tier_palette) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face="bold"),
    plot.title = element_text(size = 18, face = "bold", margin = margin(b = 10)),
    plot.subtitle = element_text(size = 12, margin = margin(b = 15))
  ) +
  labs(
    title = "Evidence Breakdown: Enzymes WITHOUT an EC Number",
    subtitle = "Area is proportional to entry count. Color represents the strongest evidence class present.",
    fill = "Highest Evidence Class"
  )

print(g_treemap_no_ec)


ggsave("plots/08_evidence_treemap_with_ec.png", g_treemap_with_ec, width = 10, height = 8, bg='white')
ggsave("plots/08_evidence_treemap_no_ec.png", g_treemap_no_ec, width = 10, height = 8, bg='white')
ggsave("plots/08_evidence_treemap.png", g_treemap, width = 10, height = 8, bg='white')


ggsave("plots/08_evidence_treemap_with_ec.pdf", g_treemap_with_ec, width = 10, height = 8, bg='white')
ggsave("plots/08_evidence_treemap_no_ec.pdf", g_treemap_no_ec, width = 10, height = 8, bg='white')
ggsave("plots/08_evidence_treemap.pdf", g_treemap, width = 10, height = 8, bg='white')


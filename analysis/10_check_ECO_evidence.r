# install.packages(c("arrow", "dplyr", "stringr", "tidyr", "UpSetR"))
library(arrow)
library(dplyr)
library(stringr)
library(tidyr)
library(UpSetR)
library(ggplot2)
format_k_m <- function(x) {  case_when(
  x >= 1e6 ~ sprintf("%.1f M", x / 1e6),
  TRUE     ~ sprintf("%.1f K", x / 1e3)
)
}

setwd("/media/benjamin/DATA/GitHub/uniprot-enzyme-finder")

# 1. Grab the full file paths from inside BOTH directories
parquet_files <- c(
  list.files("output/2026_01_v15/uniprot_enzymes_dataset", pattern = "\\.parquet$", full.names = TRUE),
  list.files("output/2026_01_v15/uniprot_non_enzymes_dataset", pattern = "\\.parquet$", full.names = TRUE)
)

# 2. Open the combined list of exact files as a single virtual dataset
ds <- open_dataset(parquet_files, format = "parquet")

# 3. Verify it worked
print(colnames(ds))

# ==============================================================================
# PART 1: Overall Fraction of Enzymes vs Non-Enzymes (Split by EC Presence)
# ==============================================================================
overall_fractions <- ds %>%
  mutate(
    has_ec = if_else(ec_completeness %in% c("Full", "Partial"), "With EC", "No EC"),
    is_enzyme = if_else(primary_class %in% c("Enzyme", "Suspected_Enzyme", "Pseudo_Enzyme"), "Enzyme", "Non-Enzyme")
  ) %>%
  group_by(has_ec, is_enzyme) %>%
  tally(name = "count") %>%
  collect() %>% # Bring the 4-row summary table into RAM
  group_by(has_ec) %>%
  mutate(
    total_in_ec_group = sum(count),
    fraction = count / total_in_ec_group,
    pct_label = paste0(round(fraction * 100, 1), "%")
  ) %>%
  arrange(has_ec, is_enzyme)

print("--- 1. OVERALL FRACTIONS ---")
print(overall_fractions)

# ==============================================================================
# PART 2: The "Clues" Profiling
# ==============================================================================
clues_profiling <- ds %>%
  mutate(
    has_ec = if_else(ec_completeness %in% c("Full", "Partial"), "With EC", "No EC"),
    is_enzyme = if_else(primary_class %in% c("Enzyme", "Suspected_Enzyme", "Pseudo_Enzyme"), "Enzyme", "Non-Enzyme")
  ) %>%
  group_by(has_ec, is_enzyme) %>%
  summarize(
    total = n(),
    count_reaction = sum(if_else(flag_has_reaction == "TRUE", 1, 0), na.rm = TRUE),
    count_active_site = sum(if_else(flag_has_active_site == "TRUE", 1, 0), na.rm = TRUE),
    count_kinetics = sum(if_else(flag_has_kinetics == "TRUE", 1, 0), na.rm = TRUE),
    count_experimental = sum(if_else(flag_has_experimental == "TRUE", 1, 0), na.rm = TRUE),
    count_any_go = sum(if_else(go_term != "", 1, 0), na.rm = TRUE),
    count_any_kw = sum(if_else(keyword != "", 1, 0), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  collect() %>%
  mutate(
    pct_reaction = (count_reaction / total) * 100,
    pct_active_site = (count_active_site / total) * 100,
    pct_kinetics = (count_kinetics / total) * 100,
    pct_experimental = (count_experimental / total) * 100,
    pct_any_go = (count_any_go / total) * 100,
    pct_any_kw = (count_any_kw / total) * 100
  )

print("--- 2. CLUES PROFILING ---")
print(clues_profiling)


# ==============================================================================
# PART 3: Isolating the "Uncharacterized" Enzymes
# ==============================================================================
uniprot_data <- ds %>%
  filter(
    #ec_completeness == "None",
    #flag_has_reaction == "FALSE",
    primary_class %in% c("Enzyme", "Suspected_Enzyme", "Pseudo_Enzyme")
  ) %>%
  mutate(
    has_active_site    = if_else(flag_has_active_site == "TRUE", TRUE, FALSE),
    has_experimental   = if_else(flag_has_experimental == "TRUE", TRUE, FALSE),
    has_catalytic_go   = if_else(str_detect(go_term, "catalytic|activity"), TRUE, FALSE),
    has_enzyme_keyword = if_else(str_detect(keyword, "ligase|kinase|hydrolase|transferase|oxidoreductase|lyase|isomerase"), TRUE, FALSE)
  ) %>%
  group_by(
    primary_class,
    has_active_site,
    has_experimental,
    has_catalytic_go,
    has_enzyme_keyword
  ) %>%
  tally(name = "count") %>%
  collect() %>%
  arrange(primary_class, desc(count))

print("--- UNCHARACTERIZED ENZYMES: EVIDENCE OVERLAP ---")
print(uncharacterized_combinations)

# ==============================================================================
# PART 4: UpSetR Plot
# ==============================================================================
# Prepare the data for UpSetR
upset_data <- uncharacterized_combinations %>%
  ungroup() %>%
  dplyr::select(-has_experimental) %>%
  mutate(across(starts_with("has_"), as.integer)) %>%
  uncount(count) %>%
  ungroup() %>%
  select(-primary_class) %>%
  as.data.frame()

# Plot the intersections
upset(
  upset_data,
  sets = c("has_active_site", "has_catalytic_go", "has_enzyme_keyword"),
  keep.order = TRUE,
  order.by = "freq",
  number.angles = 0,
  text.scale = 1.5,
  mainbar.y.label = "Number of Suspected enzymes",
  sets.x.label = "Total per Clue"
)

## EC completeness
ec_counts <- ds %>%
  group_by(ec_completeness) %>%
  tally(name = "count") %>%
  collect() # Pulls only the 3-row summary table into RAM

# Let's take a look at the summary
print(ec_counts)

ec_plot_data <- ec_counts %>%
  filter(!is.na(ec_completeness) & ec_completeness != "") %>%
  mutate(
    percentage = count / sum(count) * 100,
    label = paste0(format(round(percentage, 1), nsmall = 1), "%", "\n",format_k_m(count),""),
    ec_completeness = factor(ec_completeness, levels=c('None',"Partial","Full"), ordered = T, labels=c("No EC","Partial EC", "EC full"))
  )

# 2. Build the Donut Plot
g1 = ggplot(ec_plot_data, aes(x = 2, y = count, fill = ec_completeness)) +
  geom_bar(stat = "identity", color = "white", size = 0.5) +
  coord_polar(theta = "y", start = 0) +
  xlim(0.5, 3.5) +
  geom_text(aes(label = ec_completeness, x=3, color=ec_completeness), position = position_stack(vjust = 0.5), size = 5, fontface = "bold") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4, fontface = "bold", color='white') +
  scale_fill_brewer(palette = "Blues", aesthetics = c('color','fill')) +
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 16, face = "italic")) +
  labs(title=sprintf('N=%s Uniprot entries',format_k_m(sum(ec_counts$count))),
       subtitle="Uniprot Release 2026_01 (SwissProt + TrEMBL)")

print(g1)



cross_counts <- ds %>%
  group_by(flag_has_ec, flag_has_reaction) %>%
  tally(name = "count") %>%
  collect()

cross_table <- cross_counts %>%
  pivot_wider(
    names_from = flag_has_reaction,
    values_from = count,
    names_prefix = "Has_Reaction_"
  ) %>%
  arrange(flag_has_ec)

print(cross_table)

# Prepare the data for plotting (calculate percentages within each EC group)
plot_cross_data <- cross_counts %>%
  filter(!is.na(flag_has_ec) & flag_has_ec != "") %>%
  group_by(flag_has_ec) %>%
  mutate(
    total_in_group = sum(count),
    percentage = count / total_in_group * 100,
    # Create a clean label with the percentage
    label = paste0(format(round(percentage, 0), nsmall = 0), "%")
  )

# Build the Stacked Bar Chart
g2 = ggplot(plot_cross_data, aes(x = flag_has_ec, y = count, fill = flag_has_reaction)) +
  geom_bar(stat = "identity", position = "stack", width = 0.3) +
  # Add percentage labels inside the stacked bars
  geom_text(
    aes(label = label, color=flag_has_reaction),
    position = position_stack(vjust = 0.5),
    size = 4, vjust=-2,
    fontface = "bold",  show.legend = F
  ) +
  geom_text(
    aes(label = format_k_m(count)),
    position = position_stack(vjust = 0.5),
    size = 4, vjust=0,
    fontface = "bold",
  ) +

  theme_minimal() +
  scale_fill_brewer(palette = 4,type = 'div',direction=-1, aesthetics = c('fill','color')) +
  scale_y_continuous(labels=scales::label_number(scale_cut = scales::cut_short_scale())) +
  facet_wrap(~flag_has_ec, scales='free_y',ncol=1) +
  theme(
    panel.grid.major.y =  element_blank(), panel.grid.minor.x =  element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "top", aspect.ratio = 0.2
  ) +
  labs(
    x = "Has EC Number",
    y = "Number of Entries"
  ) +
  coord_flip(clip='off')

g2

library(patchwork)
g1 + g2


library(arrow)
library(dplyr)
library(stringr)

# 1. Define the regex patterns for the ECO codes based on UniProt classifications
# Small-scale/Targeted Experimental (EXP, IDA, IPI, IMP, IGI, IEP)
rx_exp_small <- "ECO:0000269|ECO:0000314|ECO:0000315|ECO:0000316|ECO:0000270|ECO:0000353"

# High-Throughput Experimental (HTP, HDA, HMP, HGI, HEP)
rx_exp_htp <- "ECO:0007005|ECO:0007001|ECO:0007003|ECO:0007007|ECO:0007002"

# Curated/Manual Inference (e.g., sequence similarity, ISS, IBA, phylogeny)
rx_manual_inf <- "ECO:0000250|ECO:0000318|ECO:0000319|ECO:0000320|ECO:0000321|ECO:0000266"

# Automatic / Electronic (IEA, Automatic Assertion)
rx_auto <- "ECO:0000501|ECO:0000256"

# 2. Run the Arrow Query
go_evidence_summary <- ds %>%
  filter(primary_class %in% c("Enzyme", "Suspected_Enzyme", "Pseudo_Enzyme")) %>%
  mutate(
    # Create the EC Presence Flag
    has_ec = if_else(ec_completeness %in% c("Full", "Partial"), "With EC", "No EC"),

    # Flag the presence of specific GO Evidence codes in the string
    go_small_exp = if_else(str_detect(go_term, rx_exp_small), TRUE, FALSE),
    go_htp_exp   = if_else(str_detect(go_term, rx_exp_htp), TRUE, FALSE),
    go_man_inf   = if_else(str_detect(go_term, rx_manual_inf), TRUE, FALSE),
    go_auto      = if_else(str_detect(go_term, rx_auto), TRUE, FALSE),
    has_any_go   = if_else(go_term != "", TRUE, FALSE)
  ) %>%
  mutate(
    # Hierarchically assign the HIGHEST level of GO evidence the protein possesses
    highest_go_evidence = case_when(
      go_small_exp ~ "1_Experimental_Small_Scale",
      go_htp_exp   ~ "2_Experimental_High_Throughput",
      go_man_inf   ~ "3_Manual_Inference",
      go_auto      ~ "4_Automatic_Inference",
      has_any_go   ~ "5_Other_Evidence",
      TRUE         ~ "6_No_GO_Terms"
    )
  ) %>%
  # Cross the Class, EC status, and GO Evidence
  group_by(primary_class, has_ec, highest_go_evidence) %>%
  tally(name = "protein_count") %>%
  collect() %>% # Bring the summary table into RAM
  arrange(primary_class, has_ec, highest_go_evidence) %>%
  # Calculate percentages within each Class + EC group
  group_by(primary_class, has_ec) %>%
  mutate(
    total_in_group = sum(protein_count),
    pct_of_group = round((protein_count / total_in_group) * 100, 2)
  )

# Print the resulting matrix
print("--- GO EVIDENCE CROSSED WITH EC NUMBERS ---")
print(go_evidence_summary, n=100)


# 2. Ask Arrow to group and count identical GO strings (This is blazing fast and saves RAM)
raw_go_counts <- ds %>%
  filter(
    primary_class == "Enzyme",
    ec_completeness == "None",
    go_term != ""
  ) %>%
  group_by(go_term) %>%
  tally(name = "protein_count") %>%
  collect()

# 3. Extract, unnest, and count the individual ECO codes in RAM
eco_breakdown <- raw_go_counts %>%
  mutate(
    # Extract all ECO codes from the string
    eco_codes = str_extract_all(go_term, "ECO:[0-9]{7}"),
    # Keep only UNIQUE codes per protein (so if a protein has 5 ECO:0000501 annotations, it only counts as 1 protein)
    eco_codes = lapply(eco_codes, unique)
  ) %>%
  # Unnest the list so each unique ECO code gets its own row
  unnest(eco_codes) %>%
  # Sum the protein counts for each individual ECO code
  group_by(eco_codes) %>%
  summarize(
    total_proteins_with_code = sum(protein_count),
    .groups = "drop"
  ) %>%
  arrange(desc(total_proteins_with_code))


library(arrow)
library(dplyr)
library(stringr)
library(tidyr)

# ... [Keep Steps 1, 2, and 3 exactly the same as before] ...

# ==============================================================================
# 4. Parse the Official ECO Ontology File
# ==============================================================================
# Update this path to where the eco.obo file is located inside your extracted folder
obo_file_path <- "data/eco-basic.obo"

# Read the file as raw text lines
obo_lines <- readLines(obo_file_path)

# Build the mapping table dynamically using tidyverse
eco_mapping <- tibble(line = obo_lines) %>%
  # Keep only the lines that define an ID or a Name
  filter(str_detect(line, "^id: ECO:|^name: ")) %>%
  mutate(
    # Extract the ID and Name into their own columns
    eco_codes = if_else(str_starts(line, "id: "), str_remove(line, "id: "), NA_character_),
    evidence_category = if_else(str_starts(line, "name: "), str_remove(line, "name: "), NA_character_)
  ) %>%
  # Carry the ID down to the adjacent name row
  fill(eco_codes, .direction = "down") %>%
  # Drop the "id:" rows now that the name row has the ID attached
  filter(!is.na(evidence_category)) %>%
  select(eco_codes, evidence_category) %>%
  # Ensure no duplicates (in case of typedefs or obsolete terms)
  distinct(eco_codes, .keep_all = TRUE)

# ==============================================================================
# 5. Join the official names to your dataset
# ==============================================================================
final_eco_stats <- eco_breakdown %>%
  left_join(eco_mapping, by = "eco_codes") %>%
  mutate(
    # If an ECO code somehow isn't in the OBO file, flag it
    evidence_category = replace_na(evidence_category, "Unknown Code / Obsolete"),
    pct_of_no_ec_enzymes = round((total_proteins_with_code / sum(raw_go_counts$protein_count)) * 100, 2)
  )

print("--- INDIVIDUAL ECO CODE BREAKDOWN (NO EC ENZYMES) ---")
print(final_eco_stats, n=50)




library(arrow)
library(dplyr)
library(stringr)

# 1. Filter the dataset for entries that specifically rely on ECO:0007669
automated_enzymes <- ds %>%
  filter(
    # Only look at the suspected/pseudo groups
    primary_class %in% c("Suspected_Enzyme", "Pseudo_Enzyme"),
    # Only look at proteins where the GO term contains the automatic assertion code
    str_detect(go_term, "ECO:0007669")
  ) %>%
  mutate(
    has_active_site = if_else(flag_has_active_site == "TRUE", "Has Active Site", "No Active Site")
  ) %>%
  group_by(primary_class, has_active_site) %>%
  tally(name = "protein_count") %>%
  collect() %>%
  ungroup() %>%
  mutate(
    total_automated = sum(protein_count),
    pct_of_automated_pool = round((protein_count / total_automated) * 100, 2)
  ) %>%
  arrange(primary_class, desc(protein_count))

print("--- RELIABILITY OF ECO:0007669 (AUTOMATIC ASSERTIONS) ---")
print(automated_enzymes)

# Calculate the minimum False Positive Rate
false_positives <- automated_enzymes %>%
  filter(primary_class == "Pseudo_Enzyme") %>%
  pull(protein_count) %>%
  sum()

total_pool <- automated_enzymes$total_automated[1]
fpr <- round((false_positives / total_pool) * 100, 2)

cat("\n========================================================\n")
cat("ESTIMATED FALSE POSITIVE RATE (Homology without Catalytic Function):\n")
cat("At least", fpr, "% of ECO:0007669 automated assertions are likely false positives\n")
cat("because they lack required catalytic residues (Pseudo-enzymes).\n")
cat("========================================================\n")

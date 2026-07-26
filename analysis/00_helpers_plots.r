# ==============================================================================
# 00_helpers.R
# Shared functions, palettes, themes, and dataset targets for UniProt plotting
# ==============================================================================
library(dplyr)
library(ggplot2)
library(stringr)
library(arrow)

VERSION = "v30"
RELEASE <- "2026_01"
OUTPUT_DIR = sprintf("output/%s_%s/",RELEASE,VERSION)
# Ensure the plots directory exists
PLOT_DIR = sprintf("plots/%s_%s/",RELEASE,VERSION)
dir.create(PLOT_DIR,showWarnings = F,recursive = T)

# --- 0. Uniprot proteomes ---
check_file <- function(file, organism) {
  if (!file.exists(file))
    stop("No local ", organism," UniProt reference proteome found and download failed.",call. = FALSE)
}

URL_SCER = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000002311/UP000002311_559292.fasta.gz"
URL_HSAP = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz"
FILE_SCER = "~/Downloads/UP000002311_559292.fasta.gz"
FILE_HSAP = "~/Downloads/UP000005640_9606.fasta.gz"
try(
  download.file(URL_SCER, FILE_SCER, mode = "wb"),
  download.file(URL_HSAP, FILE_HSAP, mode = "wb"),
  silent = TRUE
)

check_file(FILE_SCER, "S. cerevisiae")
check_file(FILE_HSAP, "H. sapiens")
sc_up <- Biostrings::readAAStringSet(FILE_SCER)
hs_up <- Biostrings::readAAStringSet(FILE_HSAP)

# --- 1. Dataset Targets ---
# Define the datasets we want to iterate over in every plot script
target_datasets <- list(
  "full"  = "data/plot_data/aggregated_plot_data_full.rds",
  "yeast" = "data/plot_data/aggregated_plot_data_yeast.rds"
)
# Estimate Memory Required to Load Arrow Dataset Columns
estimate_memory_gb <- function(dataset, columns, sample_size = 50000) {

  # 1. Get the total number of rows (computed instantly by Arrow)
  total_rows <- dataset %>% count() %>% pull(n)

  # 2. Lazily select the columns, take the sample, and collect into RAM
  sample_df <- dataset %>%
    select(any_of(columns)) %>%
    head(sample_size) %>%
    collect()

  # 3. Get the actual number of rows sampled (in case the dataset is small)
  actual_sample_size <- nrow(sample_df)

  if (actual_sample_size == 0) {
    warning("Sample is empty. Returning 0 GB.")
    return(0)
  }

  # 4. Measure the exact size of the sample in bytes
  size_bytes <- as.numeric(object.size(sample_df))

  # 5. Extrapolate to the entire dataset
  estimated_total_bytes <- (size_bytes / actual_sample_size) * total_rows

  # 6. Convert to Gigabytes
  estimated_gb <- round(estimated_total_bytes / (1024^3), 3)

  return(estimated_gb)
}

# --- 2. Shared Functions ---
format_k_m <- function(x) {
  case_when(
    x >= 1e6 ~ sprintf("%.1f M", x / 1e6),
    x >= 1e3 ~ sprintf("%.1f K", x / 1e3),
    TRUE ~ as.character(x)
  )
}

getmode <- function(x) {
  unique_values <- unique(x[!is.na(x)]) # Remove NA values
  unique_values[which.max(tabulate(match(x, unique_values)))]
}

# --- 3. Shared Dictionaries ---

CATALYTIC_KW = c(
  "KW-0012"= "Acyltransferase",
  "KW-0021"= "Allosteric enzyme",
  "KW-0030"= "Aminoacyl-tRNA synthetase",
  "KW-0031"= "Aminopeptidase",
  "KW-0032"= "Aminotransferase",
  "KW-0063"= "Aspartyl esterase",
  "KW-0064"= "Aspartyl protease",
  "KW-0081"= "Bacteriolytic enzyme",
  "KW-0121"= "Carboxypeptidase",
  "KW-0195"= "Cyclin (kinase)",
  "KW-0210"= "Decarboxylase",
  "KW-0223"= "Dioxygenase",
  "KW-0224"= "Dipeptidase",
  "KW-0230"= "DNA invertase",
  "KW-0239"= "DNA-directed DNA polymerase",
  "KW-0255"= "Endonuclease",
  "KW-0267"= "Excision nuclease",
  "KW-0269"= "Exonuclease",
  "KW-0326"= "Glycosidase",
  "KW-0328"= "Glycosyltransferase",
  "KW-0347"= "Helicase",
  "KW-0378"= "Hydrolase",
  "KW-0413"= "Isomerase",
  "KW-0418"= "Kinase",
  "KW-0436"= "Ligase",
  "KW-0456"= "Lyase",
  "KW-0482"= "Metalloprotease",
  "KW-0489"= "Methyltransferase",
  "KW-0503"= "Monooxygenase",
  "KW-0511"= "Multifunctional enzyme",
  "KW-0540"= "Nuclease",
  "KW-0548"= "Nucleotidyltransferase",
  "KW-0560"= "Oxidoreductase",
  "KW-0575"= "Peroxidase",
  "KW-0637"= "Prenyltransferase",
  "KW-0645"= "Protease",
  "KW-0695"= "RNA-directed DNA polymerase",
  "KW-0696"= "RNA-directed RNA polymerase",
  "KW-0697"= "Rotamase",
  "KW-0719"= "Serine esterase",
  "KW-0720"= "Serine protease",
  "KW-0723"= "Serine/threonine-protein kinase",
  "KW-0788"= "Thiol protease",
  "KW-0799"= "Topoisomerase",
  "KW-0808"= "Transferase",
  "KW-0829"= "Tyrosine-protein kinase",
  "KW-0888"= "Threonine protease",
  "KW-0904"= "Protein phosphatase",
  "KW-1278"= "Translocase"
)
evidence_abbreviations <- c(
  "GO_Catalytic"      = "GO",
  "Keyword_Catalytic" = "KW",
  "Protein_Name"      = "Name",
  "Active_Site"       = "Site",
  "Reaction_Text"     = "Rxn",
  "Rhea_ID"           = "Rhea",
  "Enzyme_DB"         = "DB",
  "Text_Mining"       = "Text",
  "EC_Number"         = "EC",
  "\\|"               = " + "
)

# --- 4. Shared Color Palettes ---
blue_palette <- c(
  "Full EC" = "#08519c",
  "Partial EC" = "#3182bd",
  "No EC" = "#c6dbef"
)

class_palette <- c(
  "Enzyme" = "#44cf32",
  "Inferred_Enzyme" = "#11bafd",
  "Suspected_Enzyme" = "#a185f8",
  "Pseudo_Enzyme" = "#df2f3f",
  "Structural_or_Transport" = "#efcd55",
  "Non_Enzyme" = "#c9c9c9"
)

tier_palette <- c(
  "1. Enzyme activity"               = "#cc1616", # Red Enzyme
  "2. High (Biochemical/Structural)" = "#2b83ba", # Strong Blue
  "3. Medium (Ontology/Homology)"    = "#abdda4", # Soft Green
  "4. Low (Lexical/Text)"            = "#fdae61", # Warning Orange
  "5. Other/Mixed"                   = "#cccccc"  # Neutral Grey
)

func_palette <- c(
  "Unclassified_Enzyme"="#DFDFDF",
  "Core_Metabolism"="#DB6D54",
  "Info_Processing"="#5D9AD1",
  "Secondary_Metabolism"="#8FC567",
  "Signaling"="#FFCF22",
  "Transport"="#efcd55",
  "Inactive_Non_Enzyme"="#efa99f"
)

# --- 5. Shared Themes & Sizes ---
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
    plot.margin = margin(3, 3, 3, 3, "mm")
  )



# Data Orders
status_order <- c("Enzyme", "Inferred_Enzyme", "Suspected_Enzyme", "Pseudo_Enzyme", "Structural_or_Transport", "Non_Enzyme")
ec_order <- c("Full EC", "Partial EC", "No EC")
tier_order <- c("Gold", "Silver", "Bronze", "Unrated")

# Colors
status_colors <- c("Enzyme" = "#5ABE6A", "Inferred_Enzyme" = "#A1DEA1",
                   "Suspected_Enzyme" = "#ECE55B", "Pseudo_Enzyme" = "#DEA645",
                   "Structural_or_Transport"="#F96F6F", "Non_Enzyme" = "#F96F6F")
ec_colors <- c("Full EC" = "#339FED", "Partial EC" = "#65CEFF", "No EC" = "#CAEFFC")
tier_colors <- c("Gold" = "#D3AF37", "Silver" = "#9BA1AB", "Bronze" = "#CD7F32", "Unrated" = "#EFCDFE")
catalytic_colors <- c("Yes" = "gray25", "No" = "gray90")
palette_text_colors <- c(status_colors, ec_colors, tier_colors, catalytic_colors)

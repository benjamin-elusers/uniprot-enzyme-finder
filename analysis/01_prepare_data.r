library(arrow)
library(dplyr)
library(stringr)
library(tidyr)
library(tidyverse)
setwd("/media/benjamin/DATA/GitHub/uniprot-enzyme-finder")

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

# 1. Connect to Dataset
parquet_files <- c(
  list.files("output/2026_01_v15/uniprot_enzymes_dataset", pattern = "\\.parquet$", full.names = TRUE),
  list.files("output/2026_01_v15/uniprot_non_enzymes_dataset", pattern = "\\.parquet$", full.names = TRUE)
)
ds <- open_dataset(parquet_files, format = "parquet")
ds_head = ds %>% head() %>% collect

ds_sc = ds %>%
  filter(organism_id %in% c('4932','559292')) %>%
  mutate(is_reference = str_detect(keyword,"Reference proteome"),
         is_enzyme = primary_class %in% c("Enzyme", "Pseudo_Enzyme", "Suspected_Enzyme")) %>%
  collect()

ds_coverage = read_tsv("output/2026_01_v15/uniprot_field_coverage.tsv")
ds_statistics = jsonlite::read_json("output/2026_01_v15/uniprot_stats_summary.json")

# 2. Base Master Summary
master_summary <- ds %>%
  mutate(
    ec_cat = case_when(
      ec_completeness == "Full" ~ "Full EC",
      ec_completeness == "Partial" ~ "Partial EC",
      TRUE ~ "No EC"
    )
  ) %>%
  group_by(ec_cat, flag_has_reaction, primary_class) %>%
  tally(name = "count") %>%
  collect()

# 3. Evidence breakdown strings for the "Dark Matter"
evidence_raw <- ds %>%
  filter(primary_class %in% c("Enzyme", "Pseudo_Enzyme", "Suspected_Enzyme")) %>%
  mutate(has_ec = if_else(ec_completeness %in% c("Full", "Partial"), "With EC", "No EC")) %>%
  group_by(has_ec, primary_class, evidence_field) %>%
  tally(name = "count") %>%
  collect()

# 4. Functional category breakdown strings for the "Dark Matter"
func_cat_raw = ds %>%
  filter(primary_class %in% c("Enzyme", "Pseudo_Enzyme", "Suspected_Enzyme")) %>%
  group_by(functional_category) %>%
  tally(name = "total") %>%
  collect()

func_cat_ec <- ds %>%
  filter(primary_class %in% c("Enzyme", "Pseudo_Enzyme", "Suspected_Enzyme")) %>%
  mutate(has_ec = if_else(ec_completeness %in% c("Full", "Partial"), "With EC", "No EC")) %>%
  group_by(has_ec, functional_category) %>%
  tally(name = "count_by_ec") %>%
  collect()

func_cat_class = ds %>%
  filter(primary_class %in% c("Enzyme", "Pseudo_Enzyme", "Suspected_Enzyme")) %>%
  mutate(has_ec = if_else(ec_completeness %in% c("Full", "Partial"), "With EC", "No EC")) %>%
  group_by(has_ec, primary_class, functional_category) %>%
  tally(name = "count_both") %>%
  collect()


ds_function = inner_join(func_cat_raw,func_cat_ec, by='functional_category') %>%
              inner_join(func_cat_class, by = c('functional_category','has_ec')) %>%
  mutate( is_enzyme = primary_class %in% c("Enzyme", "Pseudo_Enzyme", "Suspected_Enzyme")) %>%
  group_by(functional_category, has_ec, is_enzyme) %>%
  mutate(count_enzyme = sum(count_both)) %>%
  arrange(functional_category,is_enzyme,has_ec, primary_class)

# 5. Catalytic activity breakdown strings for the "Dark Matter"
cat_act_raw = ds %>%
  mutate( has_catalytic_activity = !is.na(cc_catalytic_activity) & cc_catalytic_activity != "",
          has_ec = if_else(ec_completeness %in% c("Full", "Partial"), "With EC", "No EC"),
          is_enzyme = primary_class %in% c("Enzyme", "Pseudo_Enzyme", "Suspected_Enzyme") ) %>%
  group_by(has_catalytic_activity) %>%
  tally(name = "total") %>%
  collect()

cat_act_ec = ds %>%
  mutate(has_ec = if_else(ec_completeness %in% c("Full", "Partial"), "With EC", "No EC"),
         has_catalytic_activity = !is.na(cc_catalytic_activity) & cc_catalytic_activity != "") %>%
  group_by(has_ec, has_catalytic_activity) %>%
  tally(name = "count_by_ec") %>%
  collect()

cat_act_enz = ds %>%
  mutate(has_ec = if_else(ec_completeness %in% c("Full", "Partial"), "With EC", "No EC"),
         is_enzyme = primary_class %in% c("Enzyme", "Pseudo_Enzyme", "Suspected_Enzyme"),
         has_catalytic_activity = !is.na(cc_catalytic_activity) & cc_catalytic_activity != "") %>%
  group_by(has_ec, primary_class, has_catalytic_activity) %>%
  tally(name = "count_both") %>%
  collect()

ds_catalytic = inner_join(cat_act_raw,cat_act_ec, by='has_catalytic_activity') %>%
               inner_join(cat_act_enz, by = c('has_catalytic_activity','has_ec'), relationship = "many-to-many") %>%
  mutate( is_enzyme = primary_class %in% c("Enzyme", "Pseudo_Enzyme", "Suspected_Enzyme")) %>%
  group_by(has_catalytic_activity, has_ec, is_enzyme) %>%
  mutate(count_by_enzyme = sum(count_both)) %>%
  arrange(has_catalytic_activity,is_enzyme,has_ec, primary_class)

conf = ds %>%
  filter(primary_class %in% c("Enzyme", "Pseudo_Enzyme", "Suspected_Enzyme")) %>%
  mutate(has_ec = if_else(ec_completeness %in% c("Full", "Partial"), "With EC", "No EC"),
         is_enzyme = primary_class %in% c("Enzyme", "Pseudo_Enzyme", "Suspected_Enzyme")) %>%
  group_by(has_ec, is_enzyme, confidence) %>%
  tally(name = "count") %>%
  collect()


ds_count = ds %>%
  group_by(primary_class, flag_has_reaction, ec_completeness) %>%
  count(name='n') %>%
  collect() %>%
  arrange(primary_class,ec_completeness,flag_has_reaction)

df_clues_raw <- ds %>%
  filter(primary_class %in% c("Enzyme","Pseudo_Enzyme","Suspected_Enzyme"), ec_completeness == "None") %>%
  mutate(
    cat_go       = if_else(str_detect(go_term, "catalytic|activity"), TRUE, FALSE),
    auto_eco     = if_else(str_detect(go_term, "ECO:0007669"), TRUE, FALSE),
    cat_kw       = if_else(str_detect(keyword, "ligase|kinase|hydrolase|transferase|oxidoreductase|lyase|isomerase"), TRUE, FALSE),
    act_site     = if_else(flag_has_active_site == "TRUE", TRUE, FALSE),
    exp_evidence = if_else(flag_has_experimental == "TRUE", TRUE, FALSE),
    kinetics     = if_else(flag_has_kinetics == "TRUE", TRUE, FALSE),
    enzyme_db    = if_else(xref_brenda != "" | xref_reactome != "", TRUE, FALSE),
    rhea_id      = if_else(rhea_ids != "", TRUE, FALSE)
  ) %>%
  group_by(cat_go, auto_eco, cat_kw, act_site, exp_evidence, kinetics, enzyme_db, rhea_id) %>%
  tally(name = "count") %>%
  collect()


# Get Top 5 Keywords for this subset
all_keywords <- ds %>%
  filter(primary_class %in% c("Enzyme","Pseudo_Enzyme","Suspected_Enzyme"), ec_completeness == "None", keyword != "") %>%
  select(keyword) %>%
  collect()

library(data.table)
setDT(all_keywords)
total_rows <- nrow(all_keywords)

# 1. Split and unnest keywords fast
split_dt <- all_keywords[, .(keyword = unlist(strsplit(keyword, ";", fixed = TRUE)))]

# 2. Count for all keywords
top_keywords <- split_dt[, .(count = .N), by = keyword][order(-count)]
top_keywords[, pct_of_subset := (count / total_rows) * 100]

# 3. Count for catalytic keywords only
top_catalytic <- split_dt[keyword %in% CATALYTIC_KW, .(count = .N), by = keyword][order(-count)]
top_catalytic[, pct_of_subset := (count / total_rows) * 100]


# top_keywords = all_keywords %>%
#   separate_rows(keyword, sep = ";") %>%
#   count(keyword, name = "count", sort = TRUE) %>%
#   mutate(pct_of_subset = count / nrow(all_keywords) * 100)
#
# top_catalytic = all_keywords %>%
#   separate_rows(keyword, sep = ";") %>%
#   filter(keyword %in% CATALYTIC_KW) %>%
#   count(keyword, name = "count", sort = TRUE) %>%
#   mutate(pct_of_subset = count / nrow(all_keywords) * 100)

# 4. Save to RDS
dir.create("data/plot_data", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(
  counts = ds_count,
  master = master_summary,
  yeast = ds_sc,
  evidence = evidence_raw,
  func_cat = ds_function,
  cat_act = ds_catalytic,
  top_kw = top_keywords,
  top_kw_cat = top_catalytic,
  n_total = nrow(ds)
), "data/plot_data/aggregated_plot_data.rds")
print("Aggregated data saved to data/plot_data/aggregated_plot_data.rds")

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
abbreviations <- c(
  "GO_Catalytic"      = "GO",
  "Keyword_Catalytic" = "KW",
  "Protein_Name"      = "Name",
  "Active_Site"       = "Site",
  "Reaction_Text"     = "Rxn",
  "Rhea_ID"           = "Rhea",
  "Enzyme_DB"         = "DB",
  "Text_Mining"       = "Text",
  "EC_Number"         = "EC"
)
col_key = c("accession", "id", "gene_primary","protein_name")
col_ec = c("flag_has_ec", "ec", "ec_from_name", "all_ec", "ec_family", "ec_class", "ec_completeness")
col_evi = c("confidence","keyword", "go_term","evidence_field", "evidence_value")
col_rxn = c("flag_has_reaction", "cc_catalytic_activity", "reaction_text", "substrates", "products")
col_exp = c("flag_has_kinetics","flag_has_experimental", "kinetics_km", "kinetics_vmax", "kinetics_kcat")
col_feat = c("flag_has_active_site","ft_act_site", "cc_cofactor", "ft_binding", "ft_mutagenesis","flag_is_pseudo")
col_comment = c("cc_function", "cc_pathway", "caution_activity_relevance",
                "cc_subcellular_location", "cc_subunit", "cc_interaction",
                "cc_caution", "cc_induction", "cc_similarity", "cc_disease", "cc_phenotype")
col_xref = c("rhea_ids","xref_alphafolddb","xref_brenda", "xref_reactome", "xref_string",
             "xref_ensembl","xref_kegg", "xref_unipathway")
col_org = c("organism_id", "organism_name", "kingdom", "lineage")
col_vers = c("version", "date_modified", "date_created")
ORDER_COLS = c(col_key, "is_reference", "is_enzyme",
               "primary_class", "enzyme_type","functional_category",
               col_ec, col_evi,
               col_rxn, col_exp,
               "is_passive_transport",
               "reviewed", "protein_existence", "length", "fragment",
               col_feat, col_comment, col_xref, col_org, col_vers)


# 2. Load data
dat <- readRDS("data/plot_data/aggregated_plot_data.rds")

uni_sc = dat$yeast %>% relocate(all_of(ORDER_COLS)) %>%
  mutate(ec_completeness = fct_relevel(ec_completeness,"None",after=2))


janitor::tabyl(uni_sc, primary_class, flag_has_ec, is_reference) %>%
  janitor::adorn_percentages(denominator = 'all') %>%
  janitor::adorn_pct_formatting() %>%
  janitor::adorn_ns() %>%
  knitr::kable(align='r')


uni_sc %>% filter(is_enzyme) %>%
janitor::tabyl(functional_category, is_reference) %>%
  janitor::adorn_percentages(denominator = 'col') %>%
  janitor::adorn_pct_formatting() %>%
  janitor::adorn_ns() %>%
  knitr::kable(align='r')


janitor::tabyl(uni_sc , is_enzyme, ec_completeness, is_reference)

janitor::tabyl(uni_sc,ec_completeness,is_reference)
janitor::tabyl(uni_sc,primary_class,is_reference)
janitor::tabyl(uni_sc,is_enzyme,is_reference)

list(all_yeast = uni_sc, ref_yeast_enzyme = uni_sc %>% filter(is_reference & is_enzyme), noref_yeast_enzyme = uni_sc %>% filter(!is_reference & is_enzyme)) %>%
  write_xlsx(path="output/uniprot_yeast.xlsx")

uni_sc %>% filter(!is_reference & is_enzyme & ec_family == "Non-Enzyme") %>% View()


sc_marginal <- sc_enz %>%
  #  filter(primary_class != "Suspected_Enzyme") %>%
  filter(!evidence_field %in% c("None", "")) %>%
  # If a row had count=100 and evidence="GO|EC", it becomes two rows of 100 for "GO" and "EC".
  separate_rows(evidence_field, sep = "\\|") %>%
  # Apply abbreviations to the now-individual strings
  mutate(evidence_clean = str_replace_all(evidence_field, abbreviations)) %>%
  # Group by EC status and the single evidence type, then sum the overlaps
  group_by(ec_completeness, evidence_clean) %>%
  summarize(count = n(), .groups = "drop")


sc_enz %>% filter(ec_completeness == 'None' & evidence_field != "Enzyme_DB") %>% View()


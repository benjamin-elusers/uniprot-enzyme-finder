# config/patterns.py
import re

# =============================================================================
#   REGEX PATTERNS
# =============================================================================

# Matches EC numbers in text, e.g., "(EC 1.1.1.1)"
EC_NAME_PATTERN = re.compile(r"\(EC\s+[\d\.\-n]+\)")

# Strong indicators that the protein is inactive
HARD_INACTIVE_PATTERN = re.compile(r"\b(inactive|non-catalytic|catalytic-like|pseudo)\b", re.IGNORECASE)

# Weaker indicators (could be regulatory)
SOFT_INACTIVE_PATTERN = re.compile(r"\b(suppressor|inhibitor|regulator)\b", re.IGNORECASE)

# Putative terms found in name of uncharacterized proteins
PUTATIVE_TERMS = re.compile(r"\b(putative|probable|uncharacterized|predicted)\b", re.IGNORECASE)

# "Biosynthesis of" often implies the protein IS an enzyme in that pathway
BIOSYNTH_PATTERN = re.compile(r"\bbiosynthesis of\b", re.IGNORECASE) 

# Secondary Metabolism Text Scanners
SEC_MET_PATTERNS = re.compile(
    r"\b(alkaloid|terpenoid|terpene|polyketide|antibiotic|non-ribosomal|"
    r"toxin|biosynthesis of secondary|siderophore|flavonoid|phenylpropanoid|"
    r"lignin|coumarin|cyanogenic|glucosinolate|phytoalexin|pigment|strigolactone)\b", 
    re.IGNORECASE
)

# Names that suggest Enzyme activity even without EC
ENZYME_NAME_ENDINGS = re.compile(r"(ase|synthase|kinase|transferase)$", re.IGNORECASE)


# Reaction cleaners
TRANSPORT_PATTERN = re.compile(r'\((in|out|external|internal)\)', re.IGNORECASE)
COMPARTMENT_PATTERN = re.compile(r'\[.*?\]')

# Kinetic cleaners
PUBMED_CLEAN_PATTERN = re.compile(r'\s*\(PubMed:\d+(?:,\s*PubMed:\d+)*\)', re.IGNORECASE)
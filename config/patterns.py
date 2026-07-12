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
#ENZYME_NAME_ENDINGS = re.compile(r"(ase|synthase|kinase|transferase)$", re.IGNORECASE)
# 1. Smarter Enzyme Name Scanner (Matches catalytic terms embedded anywhere)
EMBEDDED_ENZYME_NAME = re.compile(r"\b[a-zA-Z]*(?:ase|synthase|kinase|transferase)s?\b", re.IGNORECASE)
# 2. Regulatory / Accessory Subunit Penalties
REGULATORY_PREFIX = re.compile(r"\b(inhibitor|regulator|repressor|activator|receptor|scaffold|chaperone|sensor|modulator|subunit)\b", re.IGNORECASE)
# 3. Catalytic Domain Scanner (For InterPro/Pfam/SUPFAM descriptions)
CATALYTIC_DOMAINS = re.compile(r"\b(kinase|transferase|hydrolase|lyase|isomerase|ligase|oxidoreductase|polymerase|peptidase|protease)\b", re.IGNORECASE)


# Reaction cleaners
TRANSPORT_PATTERN = re.compile(r'\((in|out|external|internal)\)', re.IGNORECASE)
COMPARTMENT_PATTERN = re.compile(r'\[.*?\]')

# Kinetic cleaners
PUBMED_CLEAN_PATTERN = re.compile(r'\s*\(PubMed:\d+(?:,\s*PubMed:\d+)*\)', re.IGNORECASE)

# Kinetics Parsing Patterns
KIN_RESP_1 = re.compile(r"(?:(?:kcat|KM|Vmax)\s+(?:is|are|=)\s*)?(.*?)\s+and\s+(.*?)\s+(?:with|for|toward|against)\s+(.*?)\s+and\s+(.*?)\s+as substrates(?:,\s*respectively)?", re.IGNORECASE)
KIN_RESP_2 = re.compile(r"(?:The\s+apparent\s+)?(?:kcat|KM|Vmax)\s+of\s+(.*?)\s+and\s+(.*?)(?:\s+degradation)?\s+(?:is|are|=)\s*(.*?)\s*,?\s*and\s*(.*?)\s*,?\s*respectively", re.IGNORECASE)
KIN_MULTI_COND = re.compile(r'^(?:kcat|KM|Vmax)\s+(at\s+pH\s+[\d\.]+)\s+(?:is|are|=)\s*([~<>\s]*[\d\.]+\s*[a-zA-Z\/\-\(\)\d]+)\s*,\s*(at\s+pH\s+[\d\.]+)\s+(?:is|are|=)\s*([~<>\s]*[\d\.]+\s*[a-zA-Z\/\-\(\)\d]+)', re.IGNORECASE)
KIN_COND_EXTRACT = re.compile(r'\(\s*(in the presence of|at |pH |with |when ).+?\)$', re.IGNORECASE)

KIN_PREFIX_1 = re.compile(r'^(?:The\s+)?(?:apparent\s+)?(?:kcat|KM|Vmax)\s+(?:is|are|=|of|:)\s+', re.IGNORECASE)
KIN_PREFIX_2 = re.compile(r'^(?:The\s+)?(?:apparent\s+)?(?:kcat|KM|Vmax)\s+', re.IGNORECASE)

KIN_INV = re.compile(r'^(?:for|toward|with|against)\s+(.+?)\s+(?:is|are|=)\s*(.+)$', re.IGNORECASE)
KIN_STD = re.compile(r'^(.+?)\s+(for|toward|with|against)\s+(.+)$', re.IGNORECASE)
KIN_VAL = re.compile(r'^([~<>\s]*[\d\.,]+(?:[eE][\+\-]?\d+)?(?:x10\^\d+)?)\s*(.*)$')
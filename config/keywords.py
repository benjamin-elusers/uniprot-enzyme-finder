# config/keywords.py

# =============================================================================
#   CATALYTIC INDICATORS
# =============================================================================

# Presence of these keywords implies the protein IS an enzyme (Activity)
CATALYTIC_KW_MAP = {
    "KW-0012": "Acyltransferase",
    "KW-0021": "Allosteric enzyme", 
    "KW-0030": "Aminoacyl-tRNA synthetase",
    "KW-0031": "Aminopeptidase", 
    "KW-0032": "Aminotransferase",
    "KW-0063": "Aspartyl esterase", 
    "KW-0064": "Aspartyl protease",
    "KW-0081": "Bacteriolytic enzyme", 
    "KW-0121": "Carboxypeptidase",
    "KW-0195": "Cyclin (kinase)", 
    "KW-0210": "Decarboxylase",
    "KW-0223": "Dioxygenase", 
    "KW-0224": "Dipeptidase",
    "KW-0230": "DNA invertase", 
    "KW-0239": "DNA-directed DNA polymerase",
    "KW-0255": "Endonuclease", 
    "KW-0267": "Excision nuclease", 
    "KW-0269": "Exonuclease", 
    "KW-0326": "Glycosidase",
    "KW-0328": "Glycosyltransferase", 
    "KW-0347": "Helicase",
    "KW-0378": "Hydrolase", 
    "KW-0413": "Isomerase", 
    "KW-0418": "Kinase", 
    "KW-0436": "Ligase",
    "KW-0456": "Lyase", 
    "KW-0482": "Metalloprotease", 
    "KW-0489": "Methyltransferase", 
    "KW-0503": "Monooxygenase",
    "KW-0511": "Multifunctional enzyme", 
    "KW-0540": "Nuclease",
    "KW-0548": "Nucleotidyltransferase", 
    "KW-0560": "Oxidoreductase", 
    "KW-0575": "Peroxidase", 
    "KW-0637": "Prenyltransferase", 
    "KW-0645": "Protease", 
    "KW-0695": "RNA-directed DNA polymerase", 
    "KW-0696": "RNA-directed RNA polymerase", 
    "KW-0697": "Rotamase",
    "KW-0719": "Serine esterase", 
    "KW-0720": "Serine protease",
    "KW-0723": "Serine/threonine-protein kinase", 
    "KW-0788": "Thiol protease", 
    "KW-0799": "Topoisomerase", 
    "KW-0808": "Transferase", 
    "KW-0829": "Tyrosine-protein kinase", 
    "KW-0888": "Threonine protease", 
    "KW-0904": "Protein phosphatase", 
    "KW-1278": "Translocase"
}


# =============================================================================
#   CONTEXTUAL CATEGORIES
# =============================================================================
# Presence implies a category context, but NOT necessarily enzymatic activity
# (e.g., Toxin, Antibiotic, Metabolism)
FUNCTIONAL_CATEGORIES = {
    "Secondary_Metabolism": [
        "KW-0017",  # Alkaloid metabolism
        "KW-0044",  # Antibiotic
        "KW-0045",  # Antibiotic biosynthesis
        "KW-0078",  # Bacteriocin
        "KW-0081",  # Bacteriolytic enzyme        
        "KW-0211",  # Defensin
        "KW-0260",  # Enterotoxin
        "KW-0295",  # Fungicide
        "KW-0414",  # Isoprene biosynthesis
        "KW-0428",  # Lantibiotic
        "KW-0528",  # Neurotoxin
        "KW-0538",  # Non-ribosomal peptide synthesis
        "KW-0638",  # Phytoalexin
        "KW-0639",  # Polyketide biosynthesis
        "KW-0716",  # Siderophore
        "KW-0800",  # Toxin
        "KW-0871",  # Bacteriocin biosynthesis
        "KW-0876",  # Taxol biosynthesis
    ],
    "Signaling": [
        "KW-0021",  # Allosteric enzyme
        "KW-0071",  # Autoinducer synthesis
        "KW-0145",  # Chemotaxis
        "KW-0297",  # G-protein coupled receptor
        "KW-0342",  # GTP-binding
        "KW-0406",  # Ion transport
        "KW-0418",  # Kinase
        "KW-0497",  # Mitogen
        "KW-0527",  # Neuropeptide
        "KW-0670",  # Receptor
        "KW-0723",  # Serine/threonine-protein kinase
        "KW-0732",  # Signal
        "KW-0807",  # Transducer
        "KW-0813",  # Transport
        "KW-0829",  # Tyrosine-protein kinase
        "KW-0879",  # Wnt signaling pathway
        "KW-0904",  # Protein phosphatase
        "KW-0914",  # Notch signaling pathway
        "KW-0932",  # Cytokinin signaling pathway
    ],
    "Info_Processing": [
        "KW-0030",  # Aminoacyl-tRNA synthetase
        "KW-0156",  # Chromatin regulator
        "KW-0228",  # DNA excision
        "KW-0234",  # DNA repair
        "KW-0235",  # DNA replication
        "KW-0239",  # DNA-directed DNA polymerase
        "KW-0251",  # Elongation factor
        "KW-0255",  # Endonuclease
        "KW-0267",  # Excision nuclease
        "KW-0269",  # Exonuclease
        "KW-0347",  # Helicase
        "KW-0396",  # Initiation factor
        "KW-0436",  # Ligase
        "KW-0508",  # mRNA splicing
        "KW-0540",  # Nuclease
        "KW-0548",  # Nucleotidyltransferase
        "KW-0692",  # RNA repair
        "KW-0695",  # RNA-directed DNA polymerase
        "KW-0696",  # RNA-directed RNA polymerase
        "KW-0697",  # Rotamase
        "KW-0698",  # rRNA processing
        "KW-0799",  # Topoisomerase
        "KW-0804",  # Transcription
        "KW-0805",  # Transcription regulation
        "KW-0810",  # Translation regulation
        "KW-0819",  # tRNA processing
    ],
    "Core_Metabolism": [
        "KW-0005",  # Acetoin biosynthesis
        "KW-0028",  # Amino-acid biosynthesis
        "KW-0119",  # Carbohydrate metabolism
        "KW-0152",  # Cholesterol biosynthesis
        "KW-0275",  # Fatty acid biosynthesis
        "KW-0276",  # Fatty acid metabolism
        "KW-0310",  # Gluconeogenesis
        "KW-0324",  # Glycolysis
        "KW-0378",  # Hydrolase
        "KW-0443",  # Lipid metabolism
        "KW-0545",  # Nucleotide biosynthesis
        "KW-0546",  # Nucleotide metabolism
        "KW-0560",  # Oxidoreductase
        "KW-0595",  # Phospholipid degradation
        "KW-0602",  # Photosynthesis
        "KW-0782",  # Thiamine biosynthesis
        "KW-0808",  # Transferase
        "KW-0816",  # Tricarboxylic acid cycle
        "KW-0831",  # Vitamin metabolism
    ]
}

# =============================================================================
#   CONSTANTS
# =============================================================================

EC_CLASSES = {
    "1": "Oxidoreductases",
    "2": "Transferases",
    "3": "Hydrolases",
    "4": "Lyases",
    "5": "Isomerases",
    "6": "Ligases",
    "7": "Translocases"
}

# Keywords that are too generic to determine function alone
WEAK_KW_IDS = {
    "KW-0082": "Biosynthesis",
    "KW-0086": "Biodegradation",
    "KW-0133": "Catabolism",
    "KW-0494": "Metabolism",
    "KW-0372": "Hormone"
}

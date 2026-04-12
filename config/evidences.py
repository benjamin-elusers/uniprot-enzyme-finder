# config/evidences.py

# 1. Gold Standard: Low-Throughput Experimental
ECO_STRICT_EXP = {
    "ECO:0000269", # EXP: Inferred from Experiment
    "ECO:0000314", # IDA: Inferred from Direct Assay
    "ECO:0000353", # IPI: Inferred from Physical Interaction
    "ECO:0000315", # IMP: Inferred from Mutant Phenotype
    "ECO:0000316", # IGI: Inferred from Genetic Interaction
    "ECO:0000270", # IEP: Inferred from Expression Pattern
}

# 2. Silver Standard: High-Throughput Experimental
ECO_HTP_EXP = {
    "ECO:0006056", # HTP: High Throughput Experiment
    "ECO:0007005", # HDA: High Throughput Direct Assay
    "ECO:0007001", # HMP: High Throughput Mutant Phenotype
    "ECO:0007003", # HGI: High Throughput Genetic Interaction
    "ECO:0007007", # HEP: High Throughput Expression Pattern
}

# 3. Bronze Standard: Manual Curation & Traceable Statements
ECO_CURATED = {
    "ECO:0000304", # TAS: Traceable Author Statement
    "ECO:0000305", # IC: Inferred by Curator
    "ECO:0000250", # ISS: Inferred from Sequence/Structural Similarity
    "ECO:0000266", # ISO: Inferred from Sequence Orthology
    "ECO:0000247", # ISA: Inferred from Sequence Alignment
    "ECO:0000255", # ISM: Inferred from Sequence Model
}

# 4. Baseline: Electronic / Untraceable
ECO_AUTOMATED = {
    "ECO:0007669", # IEA: Inferred from Electronic Annotation
    "ECO:0000303", # NAS: Non-traceable Author Statement
    "ECO:0000307", # ND: No biological Data available
}
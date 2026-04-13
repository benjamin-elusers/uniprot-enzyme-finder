# config/evidences.py

# 1. Gold Standard: Low-Throughput Experimental
ECO_STRICT_EXP = {
    "ECO:0000006", # Experimental evidence (generic)
    "ECO:0000010", # protein expression evidence
    "ECO:0000023", # affinity evidence
    "ECO:0000269", # EXP: Inferred from Experiment
    "ECO:0000270", # IEP: Inferred from Expression Pattern
    "ECO:0000314", # IDA: Inferred from Direct Assay
    "ECO:0000315", # IMP: Inferred from Mutant Phenotype
    "ECO:0000316", # IGI: Inferred from Genetic Interaction
    "ECO:0000324", # imaging assay evidence
    "ECO:0000353", # IPI: Inferred from Physical Interaction
    "ECO:0001098", # microscopy evidence
    "ECO:0001115", # fluorescence evidence
    "ECO:0005504", # spectrometry evidence
}

# 2. Silver Standard: High-Throughput Experimental
ECO_HTP_EXP = {
    "ECO:0000212", # combinatorial evidence
    "ECO:0006056", # HTP: High Throughput Experiment
    "ECO:0007000", # high throughput mutant phenotypic evidence
    "ECO:0007001", # HMP: High Throughput Mutant Phenotype
    "ECO:0007003", # HGI: High Throughput Genetic Interaction
    "ECO:0007005", # HDA: High Throughput Direct Assay
    "ECO:0007007", # HEP: High Throughput Expression Pattern
    
}

# 3. Bronze Standard: Manual Curation & Traceable Statements
ECO_CURATED = {
    "ECO:0000041", # similarity evidence
    "ECO:0000361", # inferential evidence
    "ECO:0000080", # phylogenetic evidence
    "ECO:0000200", # sequence alignment evidence
    "ECO:0000201", # sequence orthology evidence
    "ECO:0000244", # Inferred from Combinatorial Evidence
    "ECO:0000247", # ISA: Inferred from Sequence Alignment
    "ECO:0000245", # RCA: Inferred from Reviewed Computational Analysis
    "ECO:0000250", # ISS: Inferred from Sequence/Structural Similarity
    "ECO:0000255", # ISM: Inferred from Sequence Model
    "ECO:0000265", # sequence orthology evidence used in automatic assertion
    "ECO:0000266", # ISO: Inferred from Sequence Orthology
    "ECO:0000304", # TAS: Traceable Author Statement
    "ECO:0000305", # IC: Inferred by Curator
    "ECO:0000312", # imported information used in manual assertion
    "ECO:0000317", # IGC: Inferred from Genomic Context
    "ECO:0000318", # IBA: Inferred from Biological aspect of Ancestor
    "ECO:0000319", # IBD: Inferred from Biological aspect of Descendant
    "ECO:0000320", # IKR: Inferred from Key Residues
    "ECO:0000321", # IRD: Inferred from Rapid Divergence
    "ECO:0001256", # imaging assay evidence used in manual assertion
    "ECO:0007744"
}

 = ECO:0000321
# 4. Baseline: Electronic / Untraceable
ECO_AUTOMATED = {
    "ECO:0000213", # EPA: Inferred from computational prediction
    "ECO:0000256", #  
    "ECO:0000259",  
    "ECO:0000303", # NAS: Non-traceable Author Statement
    "ECO:0000307", # ND: No biological Data available
    "ECO:0000313",
    "ECO:0007669", # IEA: Inferred from Electronic Annotation
    "ECO:0007672", # computational evidence
    "ECO:0007829",
    "ECO:0008006"
}


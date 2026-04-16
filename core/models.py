# core/models.py

COLUMNS = [
    # 1. Core Metadata
    "accession", "id", "protein_name", 
    "reviewed", "protein_existence", "length", "fragment",
    "version", "date_modified", "date_created",
    "organism_id", "organism_name", "kingdom", "lineage", "gene_primary",
    
    # 2. Classification & Confidence
    "primary_class", "enzyme_type", 
    "confidence_tier", "confidence_score", 
    
    # 3. EC & Functional Categories
    "ec", "ec_from_name", "all_ec",
    "ec_family", "ec_class", "ec_completeness",
    "functional_category", 
    "is_passive_transport",
    "keyword", "go_term",
    
    # 4. Reactions & Kinetics
    "rhea_ids", "reaction_text", "substrates", "products", 
    "kinetics_km", "kinetics_vmax", "kinetics_kcat",
    
    # 5. Features & Domains
    "ft_act_site", "ft_binding", "ft_mutagenesis",
    "caution_activity_relevance",
    
    # 6. Comments (CC)
    "cc_catalytic_activity", "cc_cofactor", "cc_pathway", "cc_function",
    "cc_subunit", "cc_subcellular_location", "cc_induction", "cc_similarity",
    "cc_caution", "cc_interaction", "cc_disease", "cc_phenotype",

    # 7. Raw & Parsed Evidence Maps
    "raw_evidence_map", "detailed_evidence_map", 
    "evidence_field", "evidence_value", 
    
    # 8. Precise Feature ECO Tracking
    "ec_eco_rk", "ec_eco_id",
    "catalytic_activity_eco_rk", "catalytic_activity_eco_id",
    "go_term_eco_rk", "go_term_eco_id",
    "active_site_eco_rk", "active_site_eco_id",
    "cofactor_eco_rk", "cofactor_eco_id",
    "keyword_eco_rk", "keyword_eco_id",
    "rhea_id_eco_rk", "rhea_id_eco_id",
    "kinetics_eco_rk", "kinetics_eco_id",

    # 9. Boolean Flags
    "flag_has_ec", "flag_is_multifunctional", "flag_is_promiscuous", "flag_is_pseudo",
    "flag_has_enzyme_experimental_evidence", 
    "flag_has_reaction", "flag_has_kinetics", "flag_has_active_site", 
    "flag_has_experimental", 
    
    # 10. Cross-References
    "xref_alphafolddb", "xref_brenda", "xref_reactome", 
    "xref_string", "xref_ensembl", "xref_kegg", "xref_unipathway"
]

def build_parquet_record(entry, feats, primary_class, conf, cat, fam, cls, comp, row_metadata):
    meta = entry.get('meta', {})
    text = entry.get('text', {})
    lists = entry.get('lists', {})
    ec_list = lists.get('ec_numbers', [])
    
    return {
        "accession": meta.get('accession', ''),
        "id": meta.get('id', ''),
        "protein_name": text.get('protein_name', ''),
        "reviewed": "true" if meta.get('reviewed') else "false",
        "protein_existence": meta.get('protein_existence', ''),
        "length": meta.get('length', '0'),
        "fragment": "true" if meta.get('is_fragment') else "false",
        "version": meta.get('version', ''),
        "date_modified": meta.get('date_modified', ''),
        "date_created": meta.get('date_created', ''),
        "organism_id": text.get('organism_id', ''),
        "organism_name": text.get('organism_name', ''),
        "kingdom": text.get('kingdom', ''),
        "lineage": text.get('lineage', ''),
        "gene_primary": text.get('gene_primary', ''),
        
        "primary_class": primary_class,
        "enzyme_type": primary_class, 
        "confidence_tier": row_metadata.get("confidence_tier", "Unrated"),
        "confidence_score": str(row_metadata.get("confidence_score", 0)), 
        
        "ec": ";".join(lists.get('ec_db', [])),
        "ec_from_name": ";".join(lists.get('ec_from_name', [])),
        "all_ec": ";".join(ec_list),
        "ec_family": fam, 
        "ec_class": cls, 
        "ec_completeness": comp,

        "functional_category": cat, 
        "is_passive_transport": str(row_metadata.get("flag_is_passive_transport", False)).upper(),
        "keyword": ";".join(lists.get('kw_text', lists.get('kw_ids', []))),
        "go_term": ";".join(lists.get('go_text', lists.get('go_ids', []))),
        
        "rhea_ids": ";".join(lists.get('rhea_ids', [])),
        "reaction_text": text.get('catalytic_activity', ''),
        "substrates": "||".join(lists.get('substrates', [])),
        "products": "||".join(lists.get('products', [])),
        "kinetics_km": row_metadata.get("kinetics_km") or "",
        "kinetics_vmax": row_metadata.get("kinetics_vmax") or "",
        "kinetics_kcat": row_metadata.get("kinetics_kcat") or "",
        
        "ft_act_site": ";".join(lists.get('active_sites', [])),
        "ft_binding": ";".join(lists.get('binding_sites', [])),
        "ft_mutagenesis": ";".join(lists.get('mutagenesis', [])),
        "caution_activity_relevance": "true" if feats.get('caution_activity_relevance') else "false",
        
        "cc_catalytic_activity": text.get('catalytic_activity', ''),
        "cc_cofactor": text.get('cofactor', ''),
        "cc_pathway": text.get('pathway', ''),
        "cc_function": text.get('function', ''),
        "cc_subunit": text.get('subunit', ''),
        "cc_subcellular_location": text.get('subcellular_location', ''),
        "cc_induction": text.get('induction', ''),
        "cc_similarity": text.get('similarity', ''),
        "cc_caution": text.get('caution', ''),
        "cc_interaction": text.get('interaction', ''),
        "cc_disease": text.get('disease', ''),
        "cc_phenotype": text.get('phenotype', ''),

        "raw_evidence_map": row_metadata.get("raw_evidence_map", "{}"),
        "detailed_evidence_map": row_metadata.get("detailed_evidence_map", "[]"),
        "evidence_field": row_metadata.get('evidence_field', ''),
        "evidence_value": row_metadata.get('evidence_value', ''),
        
        "ec_eco_rk": row_metadata.get("ec_eco_rk", "NONE"),
        "ec_eco_id": row_metadata.get("ec_eco_id", "None"),
        "catalytic_activity_eco_rk": row_metadata.get("catalytic_activity_eco_rk", "NONE"),
        "catalytic_activity_eco_id": row_metadata.get("catalytic_activity_eco_id", "None"),
        "go_term_eco_rk": row_metadata.get("go_term_eco_rk", "NONE"),
        "go_term_eco_id": row_metadata.get("go_term_eco_id", "None"),
        "active_site_eco_rk": row_metadata.get("active_site_eco_rk", "NONE"),
        "active_site_eco_id": row_metadata.get("active_site_eco_id", "None"),
        "cofactor_eco_rk": row_metadata.get("cofactor_eco_rk", "NONE"),
        "cofactor_eco_id": row_metadata.get("cofactor_eco_id", "None"),
        "keyword_eco_rk": row_metadata.get("keyword_eco_rk", "NONE"),
        "keyword_eco_id": row_metadata.get("keyword_eco_id", "None"),
        "rhea_id_eco_rk": row_metadata.get("rhea_id_eco_rk", "NONE"),
        "rhea_id_eco_id": row_metadata.get("rhea_id_eco_id", "None"),
        "kinetics_eco_rk": row_metadata.get("kinetics_eco_rk", "NONE"),
        "kinetics_eco_id": row_metadata.get("kinetics_eco_id", "None"),

        "flag_has_ec": str(row_metadata.get("flag_has_ec", False)).upper(),
        "flag_is_multifunctional": str(row_metadata.get("flag_is_multifunctional", False)).upper(),
        "flag_is_promiscuous": str(row_metadata.get("flag_is_promiscuous", False)).upper(),
        "flag_is_pseudo": str(row_metadata.get("flag_is_pseudo", False)).upper(),
        "flag_has_enzyme_experimental_evidence": str(row_metadata.get("flag_has_enzyme_experimental_evidence", False)).upper(),
        "flag_has_reaction": str(row_metadata.get("flag_has_reaction", False)).upper(),
        "flag_has_kinetics": str(row_metadata.get("flag_has_kinetics", False)).upper(),
        "flag_has_active_site": str(row_metadata.get("flag_has_active_site", False)).upper(),
        "flag_has_experimental": str(row_metadata.get("flag_has_experimental", False)).upper(),
        
        "xref_alphafolddb": ";".join(lists.get('xref_af', [])),
        "xref_brenda": ";".join(lists.get('xref_brenda', [])),
        "xref_reactome": ";".join(lists.get('xref_reactome', [])),
        "xref_string": ";".join(lists.get('xref_string', [])),
        "xref_ensembl": ";".join(lists.get('xref_ensembl', [])),
        "xref_kegg": ";".join(lists.get('xref_kegg', [])),
        "xref_unipathway": ";".join(lists.get('xref_unipath', []))
    }
# core/models.py

COLUMNS = [
    "accession", "id", "protein_name", 
    "reviewed", "protein_existence", "length", "fragment",
    "version", "date_modified", "date_created",
    "organism_id", "organism_name", "kingdom", "lineage", "gene_primary",
    "primary_class", "enzyme_type", 
    "confidence_tier", "confidence_score", "raw_evidence_map",
    "exp_level_catalytic_activity", "exp_level_ec", "exp_level_go",
    "ec_family", "ec_class", "ec_completeness",
    "functional_category", 
    "flag_has_ec", "flag_has_active_site", "flag_has_kinetics", 
    "flag_has_reaction", "flag_has_experimental", "flag_is_pseudo",
    "evidence_field", "evidence_value", 
    "ec", "ec_from_name", "all_ec",
    "rhea_ids", "keyword", "go_term",
    "kinetics_km", "kinetics_vmax", "kinetics_kcat",
    "reaction_text", "substrates", "products", "is_passive_transport",
    "ft_act_site", "ft_binding", "ft_mutagenesis",
    "cc_catalytic_activity", "cc_cofactor", "cc_pathway", "cc_function",
    "cc_subunit", "cc_subcellular_location", "cc_induction", "cc_similarity",
    "cc_caution", "caution_activity_relevance", "cc_interaction",
    "cc_disease", "cc_phenotype",
    "xref_alphafolddb", "xref_brenda", "xref_reactome", 
    "xref_string", "xref_ensembl", "xref_kegg", "xref_unipathway"
]


def build_parquet_record(entry, feats, primary_class, conf, cat, fam, cls, comp, row_data):
    """Converts the parsed entry and classifier results into a flat PyArrow-compatible dict."""
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
        "confidence_tier": row_data.get("confidence_tier", "Unrated"),
        "confidence_score": str(row_data.get("confidence_score", 0)),

        "ec": ";".join(lists.get('ec_db', [])),
        "ec_from_name": ";".join(lists.get('ec_from_name', [])),
        "all_ec": ";".join(ec_list),
        "ec_family": fam, "ec_class": cls, "ec_completeness": comp,
        
        "functional_category": cat, 
        "flag_has_ec": str(row_data.get("flag_has_ec", False)).upper(),
        "flag_has_active_site": str(row_data.get("flag_has_active_site", False)).upper(),
        "flag_has_kinetics": str(row_data.get("flag_has_kinetics", False)).upper(),
        "flag_has_reaction": str(row_data.get("flag_has_reaction", False)).upper(),
        "flag_has_experimental": str(row_data.get("flag_has_experimental", False)).upper(),
        "flag_is_pseudo": str(row_data.get("flag_is_pseudo", False)).upper(),
        "is_passive_transport": str(row_data.get("flag_is_passive_transport", False)).upper(),
        
        "evidence_field": row_data.get('evidence_field', ''),
        "evidence_value": row_data.get('evidence_value', ''),
        "raw_evidence_map": row_data.get('raw_evidence_map', '{}'),
        "exp_level_catalytic_activity": row_data.get("exp_level_catalytic_activity", "NONE"),
        "exp_level_ec": row_data.get("exp_level_ec", "NONE"),
        "exp_level_go": row_data.get("exp_level_go", "NONE"),

        "rhea_ids": ";".join(lists.get('rhea_ids', [])),
        "keyword": ";".join(lists.get('kw_text', lists.get('kw_ids', []))),
        "go_term": ";".join(lists.get('go_text', lists.get('go_ids', []))),
        
        "reaction_text": text.get('catalytic_activity', ''),
        "substrates": "||".join(lists.get('substrates', [])),
        "products": "||".join(lists.get('products', [])),

        "kinetics_km": row_data.get("kinetics_km"),
        "kinetics_vmax": row_data.get("kinetics_vmax"),
        "kinetics_kcat": row_data.get("kinetics_kcat"),
        
        "ft_act_site": ";".join(lists.get('active_sites', [])),
        "ft_binding": ";".join(lists.get('binding_sites', [])),
        "ft_mutagenesis": ";".join(lists.get('mutagenesis', [])),
        
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
        
        "caution_activity_relevance": "true" if feats.get('caution_activity_relevance') else "false",
        "xref_alphafolddb": ";".join(lists.get('xref_af', [])),
        "xref_brenda": ";".join(lists.get('xref_brenda', [])),
        "xref_reactome": ";".join(lists.get('xref_reactome', [])),
        "xref_string": ";".join(lists.get('xref_string', [])),
        "xref_ensembl": ";".join(lists.get('xref_ensembl', [])),
        "xref_kegg": ";".join(lists.get('xref_kegg', [])),
        "xref_unipathway": ";".join(lists.get('xref_unipath', []))
    }

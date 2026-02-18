# core/models.py
from dataclasses import dataclass, field
from typing import List, Optional

@dataclass
class EnzymeEntry:
    """
    Represents a final classified enzyme entry for CSV output.
    """
    accession: str
    id: str
    protein_name: str
    
    # Classification Results
    enzyme_type: str            # e.g., Explicit_EC, Pseudo_Enzyme
    confidence: str             # Certified, Gold, Silver, Bronze
    functional_category: str    # Secondary_Metabolism, Signaling, etc.
    is_inactive: bool           # True if text says "inactive/suppressor"
    
    # Evidence
    ec_numbers: str             # Semicolon-separated
    catalytic_activity: str     # Extracted text reaction
    active_sites: str           # Locations
    binding_sites: str          # Locations
    keywords: str               # Semicolon-separated IDs
    go_terms: str               # Semicolon-separated IDs
    
    # Metadata
    reviewed: bool
    protein_existence: str

    def to_dict(self):
        return {
            "accession": self.accession,
            "id": self.id,
            "protein_name": self.protein_name,
            "enzyme_type": self.enzyme_type,
            "confidence": self.confidence,
            "functional_category": self.functional_category,
            "is_inactive": "true" if self.is_inactive else "false",
            "ec_numbers": self.ec_numbers,
            "catalytic_activity": self.catalytic_activity,
            "active_sites": self.active_sites,
            "binding_sites": self.binding_sites,
            "keywords": self.keywords,
            "go_terms": self.go_terms,
            "reviewed": "true" if self.reviewed else "false",
            "protein_existence": self.protein_existence
        }
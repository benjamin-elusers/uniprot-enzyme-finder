# utils/obo_loader.py
import os
import sys

def load_catalytic_go_ids(obo_path, root_id="GO:0003824"):
    """
    Loads all descendant GO IDs of 'Catalytic Activity' (GO:0003824).
    Returns a set of strings.
    """
    print(f"[Init] Loading GO DAG from {obo_path}...")
    
    if not os.path.exists(obo_path):
        print(f"[Warn] OBO file not found: {obo_path}. GO validation will be skipped.")
        return set()

    try:
        from goatools.obo_parser import GODag
        go_dag = GODag(obo_path)
        
        if root_id not in go_dag:
            print(f"[Warn] Root ID {root_id} not found in OBO file.")
            return set()
            
        # Get all children plus the root itself
        valid_ids = go_dag[root_id].get_all_children()
        valid_ids.add(root_id)
        
        # Filter out obsolete terms
        final_set = {gid for gid in valid_ids if not go_dag[gid].is_obsolete}
        print(f"[Init] Loaded {len(final_set)} catalytic GO terms.")
        return final_set

    except ImportError:
        print("[Warn] 'goatools' library not installed. Install via `pip install goatools` for GO validation.")
        return set()
    except Exception as e:
        print(f"[Error] Failed to parse OBO file: {e}")
        return set()


def load_eco_definitions(obo_file_path):
    """Utility to load ECO codes and their definitions from the eco.obo file."""
    eco_defs = {}
    try:
        with open(obo_file_path, 'r', encoding='utf-8') as f:
            current_id = None
            for line in f:
                line = line.strip()
                if line.startswith("id: ECO:"):
                    current_id = line[4:]
                elif line.startswith("name: ") and current_id:
                    eco_defs[current_id] = line[6:]
                    current_id = None
    except Exception as e:
        print(f"[Warning] Could not load ECO definitions from {obo_file_path}: {e}")
    return eco_defs

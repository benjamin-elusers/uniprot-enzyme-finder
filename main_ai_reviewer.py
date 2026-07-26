import os
import json
import asyncio
import aiohttp
import pyarrow.dataset as ds
import pyarrow as pa
import pyarrow.parquet as pq
import google.generativeai as genai
from tqdm.asyncio import tqdm

# =============================================================================
# CONFIGURATION & SEMAPHORES
# =============================================================================
# Set your API keys in your environment variables
genai.configure(api_key=os.environ.get("GEMINI_API_KEY"))

# NCBI requires an email. If you have an NCBI API key, add it to speed up limits.
ENTREZ_EMAIL = "your_email@weizmann.ac.il" 
NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "") 

VERSION_DIR = "output/2026_01_v29"
DATASET_PATH = os.path.join(VERSION_DIR, "uniprot_enzymes_dataset")
OUTPUT_FILE = os.path.join(VERSION_DIR, "extrinsic_rescues.parquet")

# --- CONCURRENCY TUNING ---
# NCBI without an API key limits to ~3 req/sec. With a key: ~10 req/sec.
NCBI_SEMAPHORE = asyncio.Semaphore(5)    
GEMINI_SEMAPHORE = asyncio.Semaphore(15) 

# =============================================================================
# LLM SCHEMA & PROMPT
# =============================================================================

# Enforcing Structured Output ensures the LLM never returns chatty text
LLM_JSON_SCHEMA = {
    "type": "object",
    "properties": {
        "is_catalytic": {"type": "boolean", "description": "True if evidence points to direct catalysis."},
        "incompatible_function_found": {"type": "boolean", "description": "True if literature strictly defines it as a scaffold, transporter, chaperone, or pseudoenzyme."},
        "base_score": {
            "type": "integer", 
            "description": "Must be exactly -50, 0, 10, 20, 30, or 50 based on the biological rigor of the evidence."
        },
        "has_convergent_evidence": {
            "type": "boolean", 
            "description": "True if multiple independent PMIDs confirm the activity."
        },
        "justification": {"type": "string", "description": "A 1-2 sentence biological justification. Mandatory if score is -50."},
        "pmids_used": {"type": "array", "items": {"type": "string"}}
    },
    "required": ["is_catalytic", "incompatible_function_found", "base_score", "has_convergent_evidence", "justification", "pmids_used"]
}

model = genai.GenerativeModel(
    model_name="gemini-2.5-pro",
    generation_config={"response_mime_type": "application/json", "response_schema": LLM_JSON_SCHEMA}
)

def build_prompt(name, accession, uniprot_context, literature_xml):
    return f"""
    You are an expert computational biochemist auditing the UniProt database to annotate reference proteomes.
    Analyze the following UniProt Database Annotations and PubMed abstracts for the protein: {name} ({accession}).
    
    Your task:
    Determine if the literature or databases provide evidence that this specific protein chain directly catalyzes a chemical reaction, OR if there is conclusive evidence that it is a non-enzyme.
    
    WARNING 1 (The Pathway Clue): 
    This protein may have little to no primary annotation but might map to metabolic databases (BioCyc, UniPathway).
    Use this as a clue, but rely on the literature and domain structure to confirm it is the ACTIVE subunit.
    
    WARNING 2 (Subunit Trap):
    Many enzymes operate as multi-protein complexes.
    If the protein is a subunit of a complex, you must determine if the literature identifies THIS SPECIFIC SUBUNIT as possessing the active site or directly performing the chemistry.
    If it is merely a structural scaffold or regulatory subunit lacking an active site, set incompatible_function_found to true and assign a score of -50.
    
    WARNING 3 (Anti-Enzyme Trap): Uncharacterized proteins are NOT confirmed non-enzymes. 
    Explicitly set incompatible_function_found to true ONLY if the protein is PROVEN to have a function incompatible with catalysis (e.g., Pseudoenzyme, pure Scaffold/Adaptor, Transporter, Chaperone). Use your best judgment for other incompatible functions.
    If you confirm it is a non-enzyme, assign a score of -50. In your justification, you MUST state the exact proven function and explicitly explain why it is incompatible with enzyme activity.
    
    SCORE GUIDE (Select the highest applicable):
    -  50: In vitro purified/recombinant assay confirms direct catalysis.
    -  30: Cell-based biochemical evidence (isotope tracing, cell lysates, metabolic flux).
    -  20: In vivo genetic evidence (knockouts/CRISPR cause specific metabolic bottlenecks).
    -  10: Structural/Domain inference (identified catalytic fold, InterPro domain, or bound substrate, but no functional assay).
    -   0: Uncharacterized or Ambiguous (Review yielded no useful catalytic data, but it is not definitively proven to be a non-enzyme).
    - -50: Confirmed Non-Enzyme (Explicit proof of an incompatible function).
    
    --- UNIPROT DATABASE EXTRACT ---
    {uniprot_context}
    
    --- PUBMED ABSTRACTS ---
    {literature_xml}
    """
# ===============================================================
# ASYNC WORKER FUNCTIONS
# ===============================================================

async def fetch_abstracts(session, pmids):
    """Fetches titles and abstracts from NCBI Entrez asynchronously."""
    if not pmids: return ""
    
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "xml",
        "email": ENTREZ_EMAIL
    }
    if NCBI_API_KEY: params["api_key"] = NCBI_API_KEY
    
    async with NCBI_SEMAPHORE:
        for attempt in range(3): # Basic retry logic to prevent pipeline crashes
            try:
                async with session.get(url, params=params) as response:
                    if response.status == 200:
                        return await response.text()
            except Exception:
                await asyncio.sleep(1)
        return ""

async def review_protein(session, row):
    """Agentic worker: Fetches literature -> Queries LLM -> Calculates Score -> Returns structured data."""
    accession = row["accession"]
    name = row["protein_name"]
    raw_ev = json.loads(row["raw_evidence_map"])
    
    # Extract unique PMIDs from the evidence map
    pmids = list(set([ev.get("pubmed") for ev in raw_ev.values() if ev.get("pubmed")]))
    if not pmids:
        return None

    # Step 1: Fetch Literature
    literature_xml = await fetch_abstracts(session, pmids)
    if not literature_xml:
        return None

    # Step 2: Format what the Python pipeline already found as context
    uniprot_context = f"""
    Catalytic Domain Found: {row.get('flag_has_catalytic_domain', False)}
    PDB Structure Found: {row.get('flag_has_pdb', False)}
    Detailed Evidence Keys: {row.get('detailed_evidence_map', '[]')}
    """

    # Step 3: Query the LLM
    prompt = build_prompt(name, accession, uniprot_context, literature_xml)

    async with GEMINI_SEMAPHORE:
        try:
            response = await model.generate_content_async(prompt)
            result = json.loads(response.text)
            
            # --- APPLY EXTRINSIC SCORING LOGIC ---
            if result.get("incompatible_function_found"):
                final_score = -50
            else:
                final_score = result.get("base_score", 0)
                if result.get("has_convergent_evidence") and final_score > 0:
                    final_score += 10 # Bonus for multiple lines of evidence
            
            # Attach the processed data to join back into the R dataframe later
            return {
                "accession": accession,
                "extrinsic_score": final_score,
                "is_catalytic": result.get("is_catalytic"),
                "incompatible_function_found": result.get("incompatible_function_found"),
                "llm_justification": result.get("justification"),
                "llm_pmids": ";".join(result.get("pmids_used", []))
            }
        except Exception as e:
            print(f"LLM Error for {accession}: {e}")
            return None

# =============================================================================
# MAIN EXECUTION
# =============================================================================

async def main():
    print("Loading Parquet Dataset...")
    dataset = ds.dataset(DATASET_PATH)
    
    # Define the Twilight Zone
    print("Filtering for Twilight Zone proteins...")
    table = dataset.to_table(columns=[
        "accession", "protein_name", "confidence_score", "raw_evidence_map", 
        "flag_has_catalytic_domain", "flag_has_pdb", "detailed_evidence_map"
    ])
    df = table.to_pandas()
    
    # Define the Twilight Zone: Ambiguous score (<50) AND at least 1 PMID to review, but <= 15 (to avoid overloading context window)
    df["pmid_count"] = df["raw_evidence_map"].apply(lambda x: len([v for v in json.loads(x).values() if v.get("pubmed")]))
    twilight_df = df[(df["confidence_score"] < 50) & (df["pmid_count"] >= 1) & (df["pmid_count"] <= 15)]
    
    records = twilight_df.to_dict('records')
    print(f"Found {len(records)} ambiguous proteins to review.")

    rescues = []
    # Use aiohttp to manage hundreds of concurrent connections
    connector = aiohttp.TCPConnector(limit=50)
    async with aiohttp.ClientSession(connector=connector) as session:
        # Run all protein reviews concurrently
        tasks = [review_protein(session, row) for row in records]
        
        # tqdm wrapper provides a live progress bar for the async tasks
        for f in tqdm(asyncio.as_completed(tasks), total=len(tasks), desc="AI Review"):
            result = await f
            # Only save the result if the LLM provided a meaningful score modification
            if result and result.get("extrinsic_score", 0) != 0:
                rescues.append(result)

    if rescues:
        print(f"Successfully processed {len(rescues)} meaningful reviews! Writing to Parquet...")
        rescue_table = pa.Table.from_pylist(rescues)
        pq.write_table(rescue_table, OUTPUT_FILE)
    else:
        print("No extrinsic score modifications applied.")

if __name__ == "__main__":
    # Handles Windows-specific asyncio loops to prevent teardown errors
    if os.name == 'nt':
        asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())
    asyncio.run(main())

from Bio.SeqUtils import gc_fraction
from Bio import SeqIO
import pandas as pd
import os
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Base directory (should be the project root)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Input/output paths
OUTPUT_CSV = os.path.join(BASE_DIR, "results", "gc_analysis", "domain_gc_content.csv")

# Define domains with adjusted coordinates for our sequences
# Note: These coordinates are based on the aligned sequences
# and may differ from the original human coordinates
DOMAINS = {
    "N-terminal": (0, 279),
    "DNA-binding": (306, 876),
    "Oligomerization": (972, 1065),
    "C-terminal": (1098, 1179)
}

# Species and their file paths (using our processed sequences)
# Species IDs to look for in the sequence headers
SPECIES_IDS = {
    "Human": "Human_TP53",
    "Pig": "Pig_TP53",
    "Rhesus_Monkey": "Rhesus_TP53",
    "Rat": "Rat_TP53",
    "Zebrafish": "Zebrafish_TP53"
}

# Single input file with all sequences
INPUT_FILE = os.path.join(BASE_DIR, "data", "processed", "all_species_aligned.fasta")

def analyze_gc_content():
    """Analyze GC content across different domains of TP53 gene"""
    try:
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
        
        results = []
        
        # Process each species
        # Read sequences from single file
        sequences = list(SeqIO.parse(INPUT_FILE, "fasta"))
        
        for species, species_id in SPECIES_IDS.items():
            logging.info(f"Processing {species} sequence...")
            
            # Find the sequence for this species
            species_seq = None
            for seq in sequences:
                if seq.id == species_id:
                    species_seq = str(seq.seq)
                    break
            
            if not species_seq:
                logging.error(f"❌ {species} sequence not found in alignment")
                continue
                
            # Analyze each domain
            for domain, (start, end) in DOMAINS.items():
                try:
                    # Extract domain sequence (adjusting for 0-based indexing)
                    domain_seq = species_seq[start:end]
                    
                    if domain_seq:
                        # Calculate GC content
                        gc = gc_fraction(domain_seq) * 100
                        
                        # Store results
                        results.append({
                            "Species": species,
                            "Domain": domain,
                            "GC_Percent": gc,
                            "Length": len(domain_seq)
                        })
                        
                        logging.info(f"✅ {species} {domain}: GC% = {gc:.2f}")
                        
                except Exception as e:
                    logging.error(f"❌ Error analyzing {species} {domain}: {str(e)}")
                    continue
        
        # Create DataFrame and save results
        df = pd.DataFrame(results)
        df.to_csv(OUTPUT_CSV, index=False)
        
        logging.info(f"✅ GC analysis saved to: {OUTPUT_CSV}")
        
        # Print summary statistics
        if not df.empty:
            logging.info("\nGC Content Summary:")
            logging.info(df.groupby(["Species", "Domain"])['GC_Percent'].mean().round(2))
            
    except Exception as e:
        logging.error(f"❌ Error in GC analysis: {str(e)}")

if __name__ == "__main__":
    analyze_gc_content()

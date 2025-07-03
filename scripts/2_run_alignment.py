import os
from Bio.Align.Applications import ClustalOmegaCommandline
import os
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Base directory (should be the project root)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Input/output paths
INPUT_FILE = os.path.join(BASE_DIR, "data", "raw_sequences", "all_species_TP53.fasta")
OUTPUT_FILE = os.path.join(BASE_DIR, "data", "processed", "all_species_aligned.clustal")

def run_alignment():
    """Run multiple sequence alignment using Clustal Omega"""
    try:
        # Check if input file exists
        if not os.path.exists(INPUT_FILE):
            logging.error(f" Input file not found: {INPUT_FILE}")
            return
            
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
        
        # Create Clustal Omega command line
        cline = ClustalOmegaCommandline(
            infile=INPUT_FILE,
            outfile=OUTPUT_FILE,
            verbose=True,
            auto=True,
            force=True  # Overwrite existing file if needed
        )
        
        logging.info(f"Running alignment with command: {cline}")
        stdout, stderr = cline()
        
        if stderr:
            logging.warning(f"Clustal Omega stderr: {stderr}")
        
        logging.info(f" Alignment saved to: {OUTPUT_FILE}")
        
    except Exception as e:
        logging.error(f" Error running alignment: {str(e)}")

if __name__ == "__main__":
    run_alignment()

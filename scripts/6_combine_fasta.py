from Bio import SeqIO
import os
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Base directory (should be the project root)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Input/output paths
OUTPUT_FILE = os.path.join(BASE_DIR, "data", "raw_sequences", "all_species_TP53.fasta")

# Species and their file paths
SPECIES_FILES = {
    "Human": "TP53_HUMAN_datasets/ncbi_dataset/data/TP53.faa",
    "Pig": "TP53_Sus scrofa (pig)_datasets-4/ncbi_dataset/data/TP53.faa",
    "Zebrafish": "tp53_Danio rerio (zebrafish)_datasets-3/ncbi_dataset/data/gene.fna",
    "Rhesus_Monkey": "TP53_Macaca mulatta (Rhesus monkey)_datasets-5/ncbi_dataset/data/TP53.faa",
    "Rat": "Tp53_Rattus norvegicus (Norway rat)_datasets-2/ncbi_dataset/data/gene.fna"
}

def combine_fasta():
    """Combine FASTA sequences from all species, keeping only the longest sequence per species"""
    try:
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
        
        # Create output file
        with open(OUTPUT_FILE, "w") as outfile:
            for species, file_path in SPECIES_FILES.items():
                full_path = os.path.join(BASE_DIR, file_path)
                
                # Check if file exists
                if not os.path.exists(full_path):
                    logging.error(f" File not found: {full_path}")
                    continue
                    
                logging.info(f"Processing {species} sequences...")
                
                # Read sequences
                try:
                    # Try both FASTA and GenBank formats
                    try:
                        records = list(SeqIO.parse(full_path, "fasta"))
                    except:
                        records = list(SeqIO.parse(full_path, "genbank"))
                    
                    # If we found sequences, get the longest one
                    if records:
                        # Get the longest sequence
                        longest_record = max(records, key=lambda r: len(str(r.seq)))
                        
                        # Write to output
                        longest_record.id = f"{species}_TP53"
                        longest_record.description = f"{species} TP53 gene"
                        SeqIO.write(longest_record, outfile, "fasta")
                        
                        logging.info(f" Added {species} sequence")
                    else:
                        logging.warning(f"No sequences found in {full_path}")
                        
                except Exception as e:
                    logging.error(f" Error processing {species}: {str(e)}")
        
        logging.info(f" Combined FASTA saved to: {OUTPUT_FILE}")
        
    except Exception as e:
        logging.error(f" Error combining FASTA files: {str(e)}")

if __name__ == "__main__":
    combine_fasta()

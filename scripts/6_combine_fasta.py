from Bio import SeqIO
import os
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Base directory (should be the project root)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Output file path
OUTPUT_FILE = os.path.join(BASE_DIR, "data", "raw_sequences", "all_species_TP53.fasta")

# Species folder mapping
SPECIES_FILES = {
    "Human": "TP53_HUMAN_datasets/ncbi_dataset/data/TP53.faa",
    "Pig": "TP53_Sus scrofa (pig)_datasets-4/ncbi_dataset/data/TP53.faa",
    "Zebrafish": "tp53_Danio rerio (zebrafish)_datasets-3/ncbi_dataset/data/gene.fna",
    "Rhesus_Monkey": "TP53_Macaca mulatta (Rhesus monkey)_datasets-5/ncbi_dataset/data/TP53.faa",
    "Rat": "Tp53_Rattus norvegicus (Norway rat)_datasets-2/ncbi_dataset/data/gene.fna"
}

def combine_fasta():
    """Combine TP53 sequences from all species into a single FASTA file"""
    try:
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
        
        with open(OUTPUT_FILE, "w") as out_handle:
            for species, rel_path in SPECIES_FILES.items():
                file_path = os.path.join(BASE_DIR, rel_path)
                
                if not os.path.exists(file_path):
                    logging.error(f"❌ File not found: {file_path}")
                    continue
                    
                logging.info(f"Processing {species} sequences...")
                for record in SeqIO.parse(file_path, "fasta"):
                    if "TP53" in record.description or "p53" in record.description.lower():
                        # Update record ID and description
                        record.id = f"{species}_TP53"
                        record.description = f"{species} TP53 gene"
                        SeqIO.write(record, out_handle, "fasta")
        
        logging.info(f"✅ Combined FASTA saved to: {OUTPUT_FILE}")
        
    except Exception as e:
        logging.error(f"❌ Error combining FASTA files: {str(e)}")

if __name__ == "__main__":
    combine_fasta()

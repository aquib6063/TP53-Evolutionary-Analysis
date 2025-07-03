import os
from Bio import SeqIO

# Base directory
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DIR = os.path.join(BASE_DIR, 'data', 'raw_sequences')

# Species list
SPECIES = [
    "Human", "Pig", "Rhesus_Monkey", "Zebrafish", "Rat"
]

def download_sequences():
    """Download TP53 sequences for all species"""
    for species in SPECIES:
        # Placeholder for actual download logic
        print(f"Downloading {species} TP53 sequence...")
        # In real implementation, this would download from NCBI
        # and save to RAW_DIR

if __name__ == "__main__":
    download_sequences()

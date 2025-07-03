import os
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo
import matplotlib.pyplot as plt
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Base directory (should be the project root)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Input/output paths
ALIGNMENT_FILE = os.path.join(BASE_DIR, "data", "processed", "all_species_aligned.fasta")
DENDROGRAM_FILE = os.path.join(BASE_DIR, "results", "figures", "tp53_dendrogram.png")

def generate_dendrogram():
    """Generate a dendrogram from the aligned sequences"""
    try:
        # Check if alignment file exists
        if not os.path.exists(ALIGNMENT_FILE):
            logging.error(f"❌ Alignment file not found: {ALIGNMENT_FILE}")
            return
            
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(DENDROGRAM_FILE), exist_ok=True)
        
        # Read the alignment
        alignment = AlignIO.read(ALIGNMENT_FILE, "fasta")
        
        # Calculate distances between sequences
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(alignment)
        
        # Create the tree
        constructor = DistanceTreeConstructor()
        tree = constructor.upgma(dm)
        
        # Plot the dendrogram
        plt.figure(figsize=(15, 10))
        Phylo.draw(tree, branch_labels=lambda c: f"{c.branch_length:.2f}")
        plt.title("TP53 Gene Evolutionary Relationships")
        plt.savefig(DENDROGRAM_FILE, dpi=300, bbox_inches='tight')
        plt.close()
        
        logging.info(f"✅ Dendrogram saved to: {DENDROGRAM_FILE}")
        
    except Exception as e:
        logging.error(f"❌ Error generating dendrogram: {str(e)}")

if __name__ == "__main__":
    generate_dendrogram()

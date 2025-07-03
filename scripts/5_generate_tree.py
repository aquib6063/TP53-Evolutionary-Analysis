from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import AlignIO
import matplotlib.pyplot as plt
import os
import logging
from Bio import Phylo

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Base directory (should be the project root)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Input/output paths
ALIGNMENT_FILE = os.path.join(BASE_DIR, "data", "processed", "all_species_aligned.fasta")
OUTPUT_FILE = os.path.join(BASE_DIR, "results", "figures", "tp53_phylogeny.png")

def generate_tree():
    """Generate a phylogenetic tree from aligned sequences"""
    try:
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
        
        # Read alignment
        logging.info("Reading alignment...")
        alignment = AlignIO.read(ALIGNMENT_FILE, "fasta")
        
        # Calculate distances
        logging.info("Calculating sequence distances...")
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(alignment)
        
        # Construct tree
        logging.info("Constructing phylogenetic tree...")
        constructor = DistanceTreeConstructor(calculator)
        tree = constructor.nj(dm)
        
        # Simplify names
        for clade in tree.find_clades():
            if clade.name:
                clade.name = clade.name.split('_')[0]
        
        # Plot tree
        logging.info("Plotting tree...")
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111)
        Phylo.draw(tree, axes=ax, do_show=False)
        plt.savefig(OUTPUT_FILE, dpi=300, bbox_inches='tight')
        
        logging.info(f"✅ Phylogenetic tree saved to: {OUTPUT_FILE}")
        
        # Print distance matrix for reference
        logging.info("\nDistance matrix:")
        for i, row in enumerate(dm):
            logging.info(f"{dm.names[i]}: {list(row)}")
            
    except Exception as e:
        logging.error(f"❌ Error generating tree: {str(e)}")

if __name__ == "__main__":
    generate_tree()

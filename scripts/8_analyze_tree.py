from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import AlignIO
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import logging
from Bio import Phylo
import numpy as np

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Base directory (should be the project root)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Input/output paths
ALIGNMENT_FILE = os.path.join(BASE_DIR, "data", "processed", "all_species_aligned.fasta")
OUTPUT_DIR = os.path.join(BASE_DIR, "results", "tree_analysis")

def analyze_tree():
    """Analyze and visualize tree topology and statistics"""
    try:
        # Create output directory if it doesn't exist
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        
        # Read alignment
        logging.info("Reading alignment...")
        alignment = AlignIO.read(ALIGNMENT_FILE, "fasta")
        
        # Calculate distances
        logging.info("Calculating sequence distances...")
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(alignment)
        
        # Get species names
        species = [name.split('_')[0] for name in dm.names]
        
        # Create distance matrix DataFrame
        dist_df = pd.DataFrame(
            np.array(dm),
            index=species,
            columns=species
        )
        
        # 1. Distance matrix heatmap
        logging.info("Creating distance matrix heatmap...")
        plt.figure(figsize=(10, 8))
        sns.heatmap(dist_df, annot=True, cmap="coolwarm", fmt='.2f',
                   cbar_kws={'label': 'Sequence Identity'})
        plt.title("Pairwise Sequence Identity Matrix")
        plt.savefig(os.path.join(OUTPUT_DIR, "distance_matrix_heatmap.png"), dpi=300, bbox_inches='tight')
        
        # 2. Distance distribution
        logging.info("Creating distance distribution plot...")
        plt.figure(figsize=(10, 6))
        distances = []
        for i in range(len(dist_df)):
            for j in range(i + 1, len(dist_df)):
                distances.append(dist_df.iloc[i, j])
        
        sns.histplot(distances, kde=True, bins=20)
        plt.title("Distribution of Pairwise Sequence Distances")
        plt.xlabel("Sequence Distance")
        plt.ylabel("Frequency")
        plt.savefig(os.path.join(OUTPUT_DIR, "distance_distribution.png"), dpi=300, bbox_inches='tight')
        
        # 3. Tree statistics
        logging.info("Calculating tree statistics...")
        constructor = DistanceTreeConstructor(calculator)
        tree = constructor.nj(dm)
        
        # Calculate branch lengths
        branch_lengths = []
        for clade in tree.find_clades():
            if clade.branch_length is not None:
                branch_lengths.append(clade.branch_length)
        
        # Calculate tree statistics
        stats = {
            "Total Branch Length": sum(branch_lengths),
            "Average Branch Length": np.mean(branch_lengths),
            "Maximum Branch Length": max(branch_lengths),
            "Minimum Branch Length": min(branch_lengths),
            "Number of Internal Nodes": len(tree.get_nonterminals()),
            "Number of Tips": len(tree.get_terminals())
        }
        
        # Save statistics to CSV
        stats_df = pd.DataFrame([stats])
        stats_df.to_csv(os.path.join(OUTPUT_DIR, "tree_statistics.csv"), index=False)
        
        # 4. Branch length distribution
        logging.info("Creating branch length distribution plot...")
        plt.figure(figsize=(10, 6))
        sns.histplot(branch_lengths, kde=True, bins=20)
        plt.title("Distribution of Branch Lengths")
        plt.xlabel("Branch Length")
        plt.ylabel("Frequency")
        plt.savefig(os.path.join(OUTPUT_DIR, "branch_length_distribution.png"), dpi=300, bbox_inches='tight')
        
        # 5. Tree topology visualization with branch lengths
        logging.info("Creating detailed tree visualization...")
        plt.figure(figsize=(12, 8))
        Phylo.draw(tree, branch_labels=lambda c: f"{c.branch_length:.2f}")
        plt.title("Phylogenetic Tree with Branch Lengths")
        plt.savefig(os.path.join(OUTPUT_DIR, "detailed_tree.png"), dpi=300, bbox_inches='tight')
        
        # Print summary
        logging.info("\nTree Analysis Summary:")
        logging.info(f"Total Branch Length: {stats['Total Branch Length']:.2f}")
        logging.info(f"Average Branch Length: {stats['Average Branch Length']:.2f}")
        logging.info(f"Number of Tips: {stats['Number of Tips']}")
        logging.info(f"Number of Internal Nodes: {stats['Number of Internal Nodes']}")
        
        logging.info(f"\nResults saved to: {OUTPUT_DIR}")
        
    except Exception as e:
        logging.error(f"❌ Error analyzing tree: {str(e)}")

if __name__ == "__main__":
    analyze_tree()

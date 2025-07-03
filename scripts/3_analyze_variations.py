import os
from Bio import AlignIO
import pandas as pd

# Base directory
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROCESSED_DIR = os.path.join(BASE_DIR, 'data', 'processed')
RESULTS_DIR = os.path.join(BASE_DIR, 'results', 'variation_analysis')

def analyze_variations():
    """Analyze sequence variations across species"""
    # Load alignment
    alignment = AlignIO.read(os.path.join(PROCESSED_DIR, 'all_species_aligned.clustal'), 'clustal')
    
    # Analyze variations
    variations = []
    for i in range(len(alignment[0])):
        column = alignment[:, i]
        if len(set(column)) > 1:  # If there's variation in this position
            variations.append({
                'position': i,
                'residues': dict(zip(alignment.names, column))
            })
    
    # Save results
    df = pd.DataFrame(variations)
    df.to_csv(os.path.join(RESULTS_DIR, 'domain_variation_summary.csv'), index=False)

if __name__ == "__main__":
    analyze_variations()

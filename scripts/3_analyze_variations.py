import os
from Bio import AlignIO
import pandas as pd
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Base directory
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Input/output paths
PROCESSED_DIR = os.path.join(BASE_DIR, 'data', 'processed')
RESULTS_DIR = os.path.join(BASE_DIR, 'results', 'variation_analysis')
ALIGNMENT_FILE = os.path.join(PROCESSED_DIR, 'all_species_aligned.fasta')
OUTPUT_CSV = os.path.join(RESULTS_DIR, 'pig_specific_variations.csv')

def analyze_variations():
    """Analyze sequence variations specific to Pig TP53"""
    try:
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
        
        # Read alignment
        alignment = AlignIO.read(ALIGNMENT_FILE, "fasta")
        
        # Find Pig sequence index
        pig_index = None
        for i, record in enumerate(alignment):
            if "Pig" in record.id:
                pig_index = i
                break
        
        if pig_index is None:
            logging.error("❌ Pig sequence not found in alignment")
            return
            
        logging.info(f"Found Pig sequence at index {pig_index}")
        
        # Analyze variations
        variations = []
        for i in range(len(alignment[0])):
            pig_base = alignment[pig_index][i]
            
            # Skip gap positions
            if pig_base == "-":
                continue
                
            # Get bases from other species
            other_bases = set()
            for j, record in enumerate(alignment):
                if j != pig_index and record[i] != "-":
                    other_bases.add(record[i])
            
            # Check if Pig has a unique base
            if pig_base not in other_bases:
                variations.append({
                    "Position": i + 1,
                    "Pig_Base": pig_base,
                    "Human_Base": alignment[0][i],
                    "Other_Bases": ",".join(sorted(other_bases)),
                    "Species": [record.id for j, record in enumerate(alignment) 
                               if j != pig_index and record[i] != "-" and record[i] != pig_base]
                })
        
        # Create DataFrame
        df = pd.DataFrame(variations)
        
        # Save results
        df.to_csv(OUTPUT_CSV, index=False)
        
        logging.info(f"Found {len(df)} pig-specific variations")
        logging.info(f"✅ Results saved to: {OUTPUT_CSV}")
        
        # Print summary
        if not df.empty:
            logging.info("\nVariation summary:")
            logging.info(df["Pig_Base"].value_counts())
            
    except Exception as e:
        logging.error(f"❌ Error analyzing variations: {str(e)}")

if __name__ == "__main__":
    analyze_variations()

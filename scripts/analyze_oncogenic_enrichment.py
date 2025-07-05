import pandas as pd
from scipy.stats import fisher_exact
import logging
from pathlib import Path
import numpy as np

def setup_logging():
    """Setup logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        filename='oncogenic_enrichment.log',
        filemode='w'
    )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

def get_cosmic_data(output_file):
    """Get mock COSMIC data for testing"""
    try:
        # Read mock data
        cosmic_df = pd.read_csv('data/mock_cosmic_variants.csv')
        
        # Save to output file
        cosmic_df.to_csv(output_file, index=False)
        
        return cosmic_df
        
    except Exception as e:
        logging.error(f"Error getting mock COSMIC data: {str(e)}")
        raise

def analyze_enrichment(pig_variants_file, cosmic_data_file):
    """Analyze enrichment of pig variants in known oncogenic positions"""
    try:
        # Load data
        pig_variants = pd.read_csv(pig_variants_file)
        cosmic_data = pd.read_csv(cosmic_data_file)
        
        # Standardize column names
        pig_variants = pig_variants.rename(columns={'Position': 'position'})
        
        # Count total variants
        total_pig_variants = len(pig_variants)
        total_cosmic_variants = len(cosmic_data)
        
        # Find overlapping variants
        overlap = pd.merge(pig_variants, cosmic_data, on="position")
        overlapping_variants = len(overlap)
        
        # Calculate Fisher's exact test
        # Create contingency table
        table = [
            [overlapping_variants, total_pig_variants - overlapping_variants],
            [total_cosmic_variants, 1000000]  # Using 1M as genome background
        ]
        
        # Perform Fisher's exact test
        _, p_value = fisher_exact(table)
        
        # Calculate enrichment ratio
        enrichment_ratio = (overlapping_variants / total_pig_variants) / (total_cosmic_variants / 1000000)
        
        # Save results
        results = {
            "overlapping_variants": overlapping_variants,
            "total_pig_variants": total_pig_variants,
            "total_cosmic_variants": total_cosmic_variants,
            "enrichment_ratio": enrichment_ratio,
            "p_value": p_value
        }
        
        return results
        
    except Exception as e:
        logging.error(f"Error in enrichment analysis: {str(e)}")
        raise

def main():
    try:
        setup_logging()
        
        # Create output directory
        output_dir = Path('results/oncogenic_enrichment')
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Get mock COSMIC data
        cosmic_file = output_dir / 'cosmic_variants.csv'
        cosmic_data = get_cosmic_data(cosmic_file)
        
        # Analyze enrichment
        logging.info("Starting enrichment analysis...")
        results = analyze_enrichment(
            'results/variant_analysis/high_impact_variants.csv',
            str(cosmic_file)
        )
        
        # Save results
        results_file = output_dir / 'enrichment_results.json'
        with open(results_file, 'w') as f:
            import json
            json.dump(results, f, indent=2)
        
        # Log summary
        logging.info(f"Enrichment analysis completed")
        logging.info(f"Results saved to: {results_file}")
        
    except Exception as e:
        logging.error(f"Error in main: {str(e)}")
        raise

if __name__ == "__main__":
    main()

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import logging
import os

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Base directory
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Input/output paths
VARIANTS_FILE = os.path.join(BASE_DIR, "results", "variation_analysis", "pig_specific_variations.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "results", "variant_analysis")

def load_variants():
    """Load and preprocess variant data"""
    logging.info("Loading and analyzing variant data...")
    variants = pd.read_csv(VARIANTS_FILE)
    
    # First calculate conservation score since it's used in impact prediction
    variants['CONSERVATION'] = variants.apply(lambda row: calculate_conservation(row), axis=1)
    
    # Then add impact prediction based on nucleotide changes and conservation
    variants['IMPACT'] = variants.apply(lambda row: predict_impact(row), axis=1)
    variants['VARIANT_TYPE'] = variants.apply(lambda row: predict_variant_type(row), axis=1)
    
    # Finally add domain information
    variants['DOMAIN'] = variants.apply(lambda row: determine_domain(row), axis=1)
    
    return variants

def predict_impact(row):
    """Predict variant impact based on nucleotide change and conservation"""
    # Get the reference base (Human base)
    ref_base = row['Human_Base']
    alt_base = row['Pig_Base']
    
    # Check for stop-gain mutations (only in coding regions)
    if ref_base in ['T', 'C'] and alt_base == 'A':
        return 'HIGH'
    
    # Check for splice site mutations (near intron/exon boundaries)
    if row['Position'] % 3 == 0 or (row['Position'] + 1) % 3 == 0:
        return 'HIGH'
    
    # Check conservation
    if row['CONSERVATION'] > 3:  # More than 3 species have same base
        return 'MODERATE'
    
    return 'LOW'

def predict_variant_type(row):
    """Predict variant type based on nucleotide change"""
    ref_base = row['Human_Base']
    alt_base = row['Pig_Base']
    
    if ref_base == alt_base:
        return 'synonymous'
    
    # Check transition vs transversion
    if (ref_base in ['A', 'G'] and alt_base in ['A', 'G']) or \
       (ref_base in ['C', 'T'] and alt_base in ['C', 'T']):
        return 'transition'
    else:
        return 'transversion'

def calculate_conservation(row):
    """Calculate conservation score based on number of species with same base"""
    # Remove quotes and split the string into a list
    other_bases_str = row['Other_Bases'].strip('"[]')
    other_bases = other_bases_str.split(',')
    
    # Clean up the bases by removing any quotes
    other_bases = [base.strip().strip('"') for base in other_bases]
    
    # Count pig base if it matches any other base
    if row['Pig_Base'] in other_bases:
        return len(other_bases) + 1  # Count pig base as well
    return len(other_bases)

def determine_domain(row):
    """Determine which domain the variant falls into"""
    position = row['Position']
    
    # Define domain boundaries (approximate positions)
    if position < 10000:  # N-terminal domain
        return 'N-terminal'
    elif position < 20000:  # DNA-binding domain
        return 'DNA-binding'
    elif position < 30000:  # Oligomerization domain
        return 'Oligomerization'
    else:  # C-terminal domain
        return 'C-terminal'

def analyze_variants(variants):
    """Perform comprehensive variant analysis"""
    try:
        # Create output directory
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        
        # 1. Impact Distribution
        impact_counts = variants['IMPACT'].value_counts()
        logging.info(f"Variant impact distribution: {impact_counts.to_dict()}")
        
        # 2. High-impact variants
        high_impact = variants[variants['IMPACT'].isin(['HIGH', 'MODERATE'])]
        logging.info(f"Number of high-impact variants: {len(high_impact)}")
        
        # 3. Pathogenicity scores
        if 'POLYPHEN' in variants.columns:
            logging.info(f"Mean PolyPhen score: {variants['POLYPHEN'].mean():.2f}")
        if 'SIFT' in variants.columns:
            logging.info(f"Mean SIFT score: {variants['SIFT'].mean():.2f}")
        
        # 4. Save results
        variants.to_csv(os.path.join(OUTPUT_DIR, 'annotated_variants.csv'), index=False)
        high_impact.to_csv(os.path.join(OUTPUT_DIR, 'high_impact_variants.csv'), index=False)
        
        # 5. Generate visualizations
        generate_visualizations(variants)
        
        logging.info("Variant analysis completed successfully!")
        
    except Exception as e:
        logging.error(f"Error in variant analysis: {str(e)}")
        raise

def generate_visualizations(variants):
    """Generate analysis visualizations"""
    try:
        # Impact distribution
        plt.figure(figsize=(10, 6))
        variants['IMPACT'].value_counts().plot(kind='bar', color='skyblue')
        plt.title('Variant Impact Distribution')
        plt.xlabel('Impact Type')
        plt.ylabel('Count')
        plt.savefig(os.path.join(OUTPUT_DIR, 'impact_distribution.png'), dpi=300, bbox_inches='tight')
        
        # PolyPhen scores
        if 'POLYPHEN' in variants.columns:
            plt.figure(figsize=(10, 6))
            sns.histplot(variants['POLYPHEN'].dropna(), bins=20, kde=True, color='purple')
            plt.title('PolyPhen Pathogenicity Scores')
            plt.xlabel('Score')
            plt.ylabel('Frequency')
            plt.savefig(os.path.join(OUTPUT_DIR, 'polyphen_scores.png'), dpi=300, bbox_inches='tight')
        
        # SIFT scores
        if 'SIFT' in variants.columns:
            plt.figure(figsize=(10, 6))
            sns.histplot(variants['SIFT'].dropna(), bins=20, kde=True, color='orange')
            plt.title('SIFT Pathogenicity Scores')
            plt.xlabel('Score')
            plt.ylabel('Frequency')
            plt.savefig(os.path.join(OUTPUT_DIR, 'sift_scores.png'), dpi=300, bbox_inches='tight')
        
        # High-impact variants by domain
        if 'Domain' in variants.columns:
            plt.figure(figsize=(12, 8))
            variants.groupby(['Domain', 'IMPACT']).size().unstack().plot(kind='bar', stacked=True)
            plt.title('High-Impact Variants by Domain')
            plt.xlabel('Domain')
            plt.ylabel('Count')
            plt.savefig(os.path.join(OUTPUT_DIR, 'domain_impact.png'), dpi=300, bbox_inches='tight')
        
        plt.close('all')
        
    except Exception as e:
        logging.error(f"Error generating visualizations: {str(e)}")
        raise

def main():
    """Main function to run analysis"""
    try:
        # Load and analyze variants
        variants = load_variants()
        analyze_variants(variants)
        
    except Exception as e:
        logging.error(f"❌ Error in main analysis: {str(e)}")
        raise

if __name__ == "__main__":
    main()

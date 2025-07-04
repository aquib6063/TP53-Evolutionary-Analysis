import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy import stats
from collections import Counter
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Base directory (should be the project root)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Input/output paths
ANNOTATED_FILE = os.path.join(BASE_DIR, "results", "variation_analysis", "vep_annotated.txt")
OUTPUT_DIR = os.path.join(BASE_DIR, "results", "variant_analysis")

def analyze_variants():
    """Analyze annotated variants with comprehensive statistical analysis"""
    try:
        # Create output directory if it doesn't exist
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        
        # Read annotated variants
        variants = pd.read_csv(ANNOTATED_FILE, sep='\t', comment='#')
        
        # 1. Impact Distribution Analysis
        logging.info("Analyzing impact distribution...")
        impact_counts = variants['IMPACT'].value_counts()
        
        # Statistical test for impact distribution
        chi2, p, dof, _ = stats.chi2_contingency([impact_counts])
        logging.info(f"Chi-squared test for impact distribution: p={p:.4f}")
        
        # Plot impact distribution
        plt.figure(figsize=(10, 6))
        impact_counts.plot(kind='bar', color='skyblue')
        plt.title('Variant Impact Distribution')
        plt.xlabel('Impact Type')
        plt.ylabel('Count')
        plt.savefig(os.path.join(OUTPUT_DIR, 'impact_distribution.png'), dpi=300, bbox_inches='tight')
        
        # 2. Pathogenicity Scores Analysis
        if 'PolyPhen' in variants.columns:
            logging.info("Analyzing PolyPhen scores...")
            polyphen_scores = variants['PolyPhen'].dropna()
            
            # Statistical tests
            mean_score = polyphen_scores.mean()
            std_dev = polyphen_scores.std()
            logging.info(f"PolyPhen scores: Mean={mean_score:.2f}, StdDev={std_dev:.2f}")
            
            # Test for normal distribution
            _, p_normal = stats.normaltest(polyphen_scores)
            logging.info(f"Normality test for PolyPhen scores: p={p_normal:.4f}")
            
            # Plot distribution
            plt.figure(figsize=(10, 6))
            sns.histplot(polyphen_scores, bins=20, kde=True, color='purple')
            plt.title('PolyPhen Pathogenicity Scores')
            plt.xlabel('Score')
            plt.ylabel('Frequency')
            plt.savefig(os.path.join(OUTPUT_DIR, 'polyphen_scores.png'), dpi=300, bbox_inches='tight')
        
        if 'SIFT' in variants.columns:
            logging.info("Analyzing SIFT scores...")
            sift_scores = variants['SIFT'].dropna()
            
            # Statistical tests
            mean_score = sift_scores.mean()
            std_dev = sift_scores.std()
            logging.info(f"SIFT scores: Mean={mean_score:.2f}, StdDev={std_dev:.2f}")
            
            # Test for normal distribution
            _, p_normal = stats.normaltest(sift_scores)
            logging.info(f"Normality test for SIFT scores: p={p_normal:.4f}")
            
            # Plot distribution
            plt.figure(figsize=(10, 6))
            sns.histplot(sift_scores, bins=20, kde=True, color='orange')
            plt.title('SIFT Pathogenicity Scores')
            plt.xlabel('Score')
            plt.ylabel('Frequency')
            plt.savefig(os.path.join(OUTPUT_DIR, 'sift_scores.png'), dpi=300, bbox_inches='tight')
        
        # 3. Domain-Specific Impact Analysis
        logging.info("Analyzing domain-specific impacts...")
        domain_impact = variants.groupby(['Domain', 'IMPACT']).size().unstack().fillna(0)
        
        # Statistical test for domain distribution
        chi2_domain, p_domain, _, _ = stats.chi2_contingency(domain_impact)
        logging.info(f"Chi-squared test for domain distribution: p={p_domain:.4f}")
        
        # Plot domain-specific impact
        plt.figure(figsize=(12, 8))
        domain_impact.plot(kind='bar', stacked=True, color=['skyblue', 'lightgreen', 'lightcoral', 'gold'])
        plt.title('Domain-Specific Impact Distribution')
        plt.xlabel('Domain')
        plt.ylabel('Count')
        plt.savefig(os.path.join(OUTPUT_DIR, 'domain_impact.png'), dpi=300, bbox_inches='tight')
        
        # 4. High-Impact Variants Analysis
        logging.info("Analyzing high-impact variants...")
        high_impact = variants[variants['IMPACT'].isin(['HIGH', 'MODERATE'])]
        
        # Statistical analysis of high-impact variants
        if len(high_impact) > 0:
            # Domain distribution of high-impact variants
            high_impact_domain = high_impact['Domain'].value_counts(normalize=True)
            logging.info(f"Domain distribution of high-impact variants: {high_impact_domain.to_dict()}")
            
            # Conservation analysis
            if 'CONSERVATION' in high_impact.columns:
                cons_scores = high_impact['CONSERVATION'].dropna()
                mean_cons = cons_scores.mean()
                logging.info(f"Mean conservation score of high-impact variants: {mean_cons:.2f}")
                
                # Plot conservation scores
                plt.figure(figsize=(10, 6))
                sns.boxplot(x='Domain', y='CONSERVATION', data=high_impact)
                plt.title('Conservation Scores of High-Impact Variants')
                plt.xlabel('Domain')
                plt.ylabel('Conservation Score')
                plt.savefig(os.path.join(OUTPUT_DIR, 'high_impact_conservation.png'), dpi=300, bbox_inches='tight')
        
        # 5. Variant Type Analysis
        logging.info("Analyzing variant types...")
        if 'VARIANT_TYPE' in variants.columns:
            variant_types = variants['VARIANT_TYPE'].value_counts()
            logging.info(f"Variant type distribution: {variant_types.to_dict()}")
            
            # Plot variant types
            plt.figure(figsize=(10, 6))
            variant_types.plot(kind='pie', autopct='%1.1f%%')
            plt.title('Variant Type Distribution')
            plt.ylabel('')
            plt.savefig(os.path.join(OUTPUT_DIR, 'variant_types.png'), dpi=300, bbox_inches='tight')
        
        # 6. Correlation Analysis
        logging.info("Analyzing correlations...")
        if 'PolyPhen' in variants.columns and 'SIFT' in variants.columns:
            # Calculate correlation
            correlation = variants[['PolyPhen', 'SIFT']].corr().iloc[0,1]
            logging.info(f"Correlation between PolyPhen and SIFT scores: {correlation:.2f}")
            
            # Plot correlation
            plt.figure(figsize=(10, 6))
            sns.scatterplot(x='PolyPhen', y='SIFT', data=variants)
            plt.title('Correlation between PolyPhen and SIFT Scores')
            plt.xlabel('PolyPhen Score')
            plt.ylabel('SIFT Score')
            plt.savefig(os.path.join(OUTPUT_DIR, 'score_correlation.png'), dpi=300, bbox_inches='tight')
        
        # Save comprehensive results
        high_impact.to_csv(os.path.join(OUTPUT_DIR, 'high_impact_variants.csv'), index=False)
        variants.to_csv(os.path.join(OUTPUT_DIR, 'annotated_variants_summary.csv'), index=False)
        
        logging.info(f"Variant analysis completed. Results saved to: {OUTPUT_DIR}")
        
    except Exception as e:
        logging.error(f"❌ Error analyzing variants: {str(e)}")

if __name__ == "__main__":
    analyze_variants()

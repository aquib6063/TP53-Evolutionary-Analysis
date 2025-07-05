import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import json
import logging

# Setup logging
def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def load_results():
    """Load enrichment analysis results"""
    try:
        with open('results/oncogenic_enrichment/enrichment_results.json', 'r') as f:
            results = json.load(f)
        
        # Load pig variants
        pig_variants = pd.read_csv('results/variant_analysis/high_impact_variants.csv')
        
        # Load cosmic variants
        cosmic_variants = pd.read_csv('results/oncogenic_enrichment/cosmic_variants.csv')
        
        return results, pig_variants, cosmic_variants
        
    except Exception as e:
        logging.error(f"Error loading results: {str(e)}")
        raise

def plot_enrichment(results, pig_variants, cosmic_variants):
    """Create enrichment visualization"""
    try:
        # Create output directory
        output_dir = Path('results/oncogenic_enrichment/figures')
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Standardize column names and data types
        cosmic_variants = cosmic_variants.rename(columns={'position': 'Position'})
        cosmic_variants['Position'] = cosmic_variants['Position'].astype(int)
        pig_variants['Position'] = pig_variants['Position'].astype(int)
        
        # Plot 1: Variant distribution across domains
        plt.figure(figsize=(10, 6))
        domain_counts = pig_variants['DOMAIN'].value_counts()
        cosmic_domain_counts = cosmic_variants['Position'].apply(
            lambda x: pig_variants[pig_variants['Position'] == x]['DOMAIN'].iloc[0]
            if x in pig_variants['Position'].values else 'Other'
        ).value_counts()
        
        # Normalize counts
        domain_counts = domain_counts / domain_counts.sum()
        cosmic_domain_counts = cosmic_domain_counts / cosmic_domain_counts.sum()
        
        # Create bar plot
        domains = list(domain_counts.index)
        pig_vals = domain_counts.values
        cosmic_vals = cosmic_domain_counts.reindex(domains, fill_value=0).values
        
        x = range(len(domains))
        width = 0.35
        
        plt.bar(x, pig_vals, width, label='Pig Variants')
        plt.bar([i + width for i in x], cosmic_vals, width, label='COSMIC Variants')
        
        plt.xlabel('Domain')
        plt.ylabel('Proportion of Variants')
        plt.title('Variant Distribution Across TP53 Domains')
        plt.xticks([i + width/2 for i in x], domains, rotation=45)
        plt.legend()
        plt.tight_layout()
        
        plt.savefig(output_dir / 'domain_distribution.png', dpi=300)
        plt.close()
        
        # Plot 2: Position-based enrichment
        plt.figure(figsize=(12, 6))
        
        # Get overlapping positions
        overlap = pd.merge(pig_variants, cosmic_variants, on='Position')
        
        # Create scatter plot
        plt.scatter(
            pig_variants['Position'],
            [0] * len(pig_variants),
            color='blue',
            alpha=0.3,
            label='Pig Variants'
        )
        plt.scatter(
            cosmic_variants['Position'],
            [0.1] * len(cosmic_variants),
            color='red',
            alpha=0.5,
            label='COSMIC Variants'
        )
        plt.scatter(
            overlap['Position'],
            [0.2] * len(overlap),
            color='purple',
            label='Overlapping Variants'
        )
        
        plt.yticks([0, 0.1, 0.2], ['Pig', 'COSMIC', 'Overlap'])
        plt.xlabel('Position')
        plt.title('Variant Position Distribution')
        plt.legend()
        plt.tight_layout()
        
        plt.savefig(output_dir / 'position_distribution.png', dpi=300)
        plt.close()
        
        # Plot 3: Enrichment ratio by domain
        plt.figure(figsize=(10, 6))
        
        # Calculate enrichment by domain
        domain_enrichment = []
        for domain in domains:
            pig_domain = pig_variants[pig_variants['DOMAIN'] == domain]
            
            # Get overlapping positions
            overlap = pd.merge(pig_domain, cosmic_variants, on='Position')
            
            if len(pig_domain) > 0:
                enrichment = (len(overlap) / len(pig_domain)) / (len(cosmic_variants) / len(cosmic_variants))
                domain_enrichment.append({
                    'domain': domain,
                    'enrichment': enrichment
                })
        
        domain_df = pd.DataFrame(domain_enrichment)
        
        # Create bar plot
        sns.barplot(data=domain_df, x='domain', y='enrichment', palette='viridis')
        plt.xlabel('Domain')
        plt.ylabel('Enrichment Ratio')
        plt.title('Enrichment of Oncogenic Variants by Domain')
        plt.xticks(rotation=45)
        plt.tight_layout()
        
        plt.savefig(output_dir / 'domain_enrichment.png', dpi=300)
        plt.close()
        
        logging.info("Visualization completed successfully")
        
    except Exception as e:
        logging.error(f"Error creating visualization: {str(e)}")
        raise

def main():
    try:
        setup_logging()
        logging.info("Starting visualization...")
        
        results, pig_variants, cosmic_variants = load_results()
        plot_enrichment(results, pig_variants, cosmic_variants)
        
        logging.info("Visualization completed successfully")
        
    except Exception as e:
        logging.error(f"Error in main: {str(e)}")
        raise

if __name__ == "__main__":
    main()

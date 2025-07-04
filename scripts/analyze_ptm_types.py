import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging

# Configure logging
logging.basicConfig(
    filename='ptm_type_analysis.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def load_ptm_data():
    """Load PTM disruption data"""
    base_dir = Path(__file__).parent.parent
    ptm_dir = base_dir / "results" / "ptm_analysis"
    
    # Load combined results
    combined_df = pd.read_csv(ptm_dir / "ptm_disruptions_combined.csv")
    
    # Load domain-specific results
    domain_dfs = []
    for file in ptm_dir.glob("ptm_disruptions_domain_*.csv"):
        df = pd.read_csv(file)
        domain_dfs.append(df)
    
    return combined_df, domain_dfs

def analyze_ptm_types(df):
    """Analyze PTM type distributions"""
    # Count PTM types by species
    human_ptms = df['human_ptm_type'].value_counts()
    pig_ptms = df['pig_ptm_type'].value_counts()
    
    # Calculate PTM type changes
    ptm_changes = df.groupby(['human_ptm_type', 'pig_ptm_type']).size()
    
    return human_ptms, pig_ptms, ptm_changes

def plot_ptm_type_distribution(human_ptms, pig_ptms, output_dir):
    """Plot PTM type distribution"""
    plt.figure(figsize=(10, 6))
    
    # Plot human PTM types
    plt.subplot(1, 2, 1)
    human_ptms.plot(kind='bar', color='blue', alpha=0.7)
    plt.title('Human PTM Types')
    plt.xlabel('PTM Type')
    plt.ylabel('Count')
    plt.xticks(rotation=45)
    
    # Plot pig PTM types
    plt.subplot(1, 2, 2)
    pig_ptms.plot(kind='bar', color='orange', alpha=0.7)
    plt.title('Pig PTM Types')
    plt.xlabel('PTM Type')
    plt.ylabel('Count')
    plt.xticks(rotation=45)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'ptm_type_distribution.png', dpi=300)
    plt.close()

def plot_ptm_changes(ptm_changes, output_dir):
    """Plot PTM type changes"""
    plt.figure(figsize=(12, 8))
    
    # Create heatmap
    ptm_changes = ptm_changes.unstack().fillna(0)
    sns.heatmap(ptm_changes, annot=True, cmap='coolwarm', fmt='g')
    plt.title('PTM Type Changes between Human and Pig')
    plt.xlabel('Pig PTM Type')
    plt.ylabel('Human PTM Type')
    
    plt.savefig(output_dir / 'ptm_type_changes.png', dpi=300)
    plt.close()

def analyze_domain_patterns(df, output_dir):
    """Analyze PTM patterns across domains"""
    # Group by domain
    domain_stats = df.groupby(['domain_start', 'domain_end']).agg({
        'human_ptm_type': pd.Series.mode,
        'pig_ptm_type': pd.Series.mode,
        'sequence_similarity': 'mean',
        'disruption_type': pd.Series.mode
    }).reset_index()
    
    # Plot similarity vs disruptions
    plt.figure(figsize=(12, 6))
    sns.scatterplot(data=domain_stats, 
                   x='sequence_similarity', 
                   y='disruption_type',
                   hue='domain_start',
                   size='domain_end')
    plt.title('Domain-wise PTM Disruption Analysis')
    plt.xlabel('Sequence Similarity')
    plt.ylabel('Disruption Type')
    plt.savefig(output_dir / 'domain_analysis.png', dpi=300)
    plt.close()
    
    return domain_stats

def main():
    try:
        # Create output directory
        output_dir = Path(__file__).parent.parent / "results" / "ptm_analysis" / "visualizations"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Load data
        combined_df, domain_dfs = load_ptm_data()
        
        # Analyze PTM types
        human_ptms, pig_ptms, ptm_changes = analyze_ptm_types(combined_df)
        
        # Generate plots
        plot_ptm_type_distribution(human_ptms, pig_ptms, output_dir)
        plot_ptm_changes(ptm_changes, output_dir)
        domain_stats = analyze_domain_patterns(combined_df, output_dir)
        
        # Save detailed analysis
        domain_stats.to_csv(output_dir / 'domain_analysis.csv', index=False)
        
        logging.info("PTM type analysis completed successfully")
        
    except Exception as e:
        logging.error(f"Error in PTM type analysis: {str(e)}")
        raise

if __name__ == "__main__":
    main()

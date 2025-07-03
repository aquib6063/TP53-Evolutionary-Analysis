import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Base directory
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VARIATION_DIR = os.path.join(BASE_DIR, 'results', 'variation_analysis')
GC_DIR = os.path.join(BASE_DIR, 'results', 'gc_content')
FIGURES_DIR = os.path.join(BASE_DIR, 'results', 'figures')

def generate_figures():
    """Generate analysis figures"""
    # Create figures directory if it doesn't exist
    os.makedirs(FIGURES_DIR, exist_ok=True)
    
    # Read data
    variations = pd.read_csv(os.path.join(VARIATION_DIR, 'domain_variation_summary.csv'))
    gc_content = pd.read_csv(os.path.join(GC_DIR, 'gc_by_domain.csv'))
    
    # Generate variation plot
    plt.figure(figsize=(10, 6))
    sns.countplot(x='position', data=variations)
    plt.title('Variation Distribution Across Sequence')
    plt.savefig(os.path.join(FIGURES_DIR, 'domain_variation_plot.png'))
    
    # Generate GC content plot
    plt.figure(figsize=(10, 6))
    sns.barplot(x='species', y='gc_content', data=gc_content)
    plt.title('GC Content by Species')
    plt.savefig(os.path.join(FIGURES_DIR, 'gc_content_changes.png'))

if __name__ == "__main__":
    generate_figures()

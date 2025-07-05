import pandas as pd
import numpy as np
import logging
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FuncFormatter
import matplotlib.patches as patches

# Setup logging
def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def create_domain_comparison_plot():
    """Create domain comparison plot using matplotlib"""
    try:
        # Create output directory
        output_dir = Path('results/interactive_plots')
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Domain data
        domains = ['N-terminal', 'DNA-binding', 'Oligomerization', 'C-terminal']
        human_gc = [55, 54, 53, 52]
        pig_gc = [55, 54, 53, 37.2]  # Note: C-terminal is significantly lower
        conservation = [89, 98, 40, 9]
        
        # Create figure
        fig, ax1 = plt.subplots(figsize=(10, 6))
        ax2 = ax1.twinx()
        
        # Create bar plots
        bar_width = 0.35
        x = np.arange(len(domains))
        
        bars1 = ax1.bar(x - bar_width/2, human_gc, bar_width, label='Human GC%', color='blue', alpha=0.7)
        bars2 = ax1.bar(x + bar_width/2, pig_gc, bar_width, label='Pig GC%', color='orange', alpha=0.7)
        
        # Create line plot for conservation
        ax2.plot(domains, conservation, 'r-', marker='o', label='Conservation %', linewidth=2)
        
        # Format axes
        ax1.set_xlabel('Protein Domain')
        ax1.set_ylabel('GC Content (%)')
        ax2.set_ylabel('Residue Conservation (%)')
        
        # Add title and grid
        plt.title('Domain-wise Comparison of GC Content and Conservation')
        plt.grid(True, alpha=0.3)
        
        # Add legend
        lines, labels = ax2.get_legend_handles_labels()
        bars, bar_labels = ax1.get_legend_handles_labels()
        ax2.legend(lines + bars, labels + bar_labels, loc='upper right')
        
        # Save PNG
        png_file = output_dir / 'domain_comparison.png'
        plt.tight_layout()
        plt.savefig(png_file, dpi=300, bbox_inches='tight')
        plt.close()
        logging.info(f"PNG saved to: {png_file}")
        
    except Exception as e:
        logging.error(f"Error creating domain comparison plot: {str(e)}")
        raise

def create_variant_position_plot():
    """Create variant position plot using matplotlib"""
    try:
        # Create output directory
        output_dir = Path('results/interactive_plots')
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Load variant data
        variants = pd.read_csv('results/variant_analysis/high_impact_variants.csv')
        
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 6))
        
        # Create scatter plot with different colors for impact levels
        color_map = {
            'HIGH': 'red',
            'MODERATE': 'orange',
            'LOW': 'green'
        }
        
        for impact, color in color_map.items():
            mask = variants['IMPACT'] == impact
            ax.scatter(
                variants[mask]['Position'],
                [1] * mask.sum(),
                color=color,
                label=f'{impact} Impact',
                alpha=0.7
            )
        
        # Add domain boundaries
        domain_boundaries = {
            'N-terminal': (1, 42),
            'DNA-binding': (94, 289),
            'Oligomerization': (323, 355),
            'C-terminal': (363, 393)
        }
        
        # Add domain rectangles and labels
        for name, (start, end) in domain_boundaries.items():
            ax.add_patch(
                patches.Rectangle(
                    (start, 0),
                    end - start,
                    2,
                    alpha=0.3,
                    color='gray'
                )
            )
            ax.text(
                (start + end) / 2,
                1.5,
                name,
                ha='center',
                fontsize=12
            )
        
        # Format axes and labels
        ax.set_xlabel('Amino Acid Position')
        ax.set_yticks([0, 1, 2])
        ax.set_yticklabels(['Low', 'Medium', 'High'])
        ax.set_ylim(-0.5, 2.5)
        
        # Add title and legend
        plt.title('Variant Positions Across TP53 Domains')
        plt.legend()
        
        # Save PNG
        png_file = output_dir / 'variant_positions.png'
        plt.tight_layout()
        plt.savefig(png_file, dpi=300, bbox_inches='tight')
        plt.close()
        logging.info(f"PNG saved to: {png_file}")
        
    except Exception as e:
        logging.error(f"Error creating variant position plot: {str(e)}")
        raise

def create_dnds_plot():
    """Create dN/dS plot using matplotlib"""
    try:
        # Create output directory
        output_dir = Path('results/interactive_plots')
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Domain data
        domains = ['N-terminal', 'DNA-binding', 'Oligomerization', 'C-terminal']
        dnds = [1.2, 0.8, 1.5, 8.49]
        
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Create bar plot
        bars = ax.bar(domains, dnds, color='skyblue', alpha=0.7)
        
        # Add neutral evolution line
        ax.axhline(y=1, color='red', linestyle='--', label='Neutral Evolution')
        
        # Add text annotation for neutral evolution
        ax.text(
            len(domains) - 0.5,
            1.05,
            'Neutral Evolution',
            color='red',
            ha='right'
        )
        
        # Format axes and labels
        ax.set_xlabel('Protein Domain')
        ax.set_ylabel('dN/dS Ratio')
        ax.set_title('dN/dS Ratios Across TP53 Domains')
        
        # Add grid
        ax.grid(True, alpha=0.3)
        
        # Save PNG
        png_file = output_dir / 'dnds_ratios.png'
        plt.tight_layout()
        plt.savefig(png_file, dpi=300, bbox_inches='tight')
        plt.close()
        logging.info(f"PNG saved to: {png_file}")
        
    except Exception as e:
        logging.error(f"Error creating dN/dS plot: {str(e)}")
        raise

def main():
    try:
        setup_logging()
        logging.info("Creating interactive plots...")
        
        create_domain_comparison_plot()
        create_variant_position_plot()
        create_dnds_plot()
        
        logging.info("Interactive plots created successfully")
        
    except Exception as e:
        logging.error(f"Error in main: {str(e)}")
        raise

if __name__ == "__main__":
    main()

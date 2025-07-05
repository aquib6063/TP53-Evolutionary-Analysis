import matplotlib.pyplot as plt
import seaborn as sns
import logging
from pathlib import Path

# Setup logging
def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def plot_selection_conservation():
    """Plot dN/dS ratios and conservation scores across TP53 domains"""
    try:
        # Create output directory
        output_dir = Path('results/figures')
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Data
        domains = ['N-terminal', 'DNA-binding', 'Oligomerization', 'C-terminal']
        dnds = [1.2, 0.8, 1.5, 8.49]
        conservation = [89, 98, 40, 9]
        
        # Plot
        fig, ax1 = plt.subplots(figsize=(10,6))
        ax2 = ax1.twinx()
        
        # Bar plot for dN/dS
        sns.barplot(x=domains, y=dnds, ax=ax1, color='skyblue', alpha=0.7)
        
        # Line plot for conservation
        sns.lineplot(x=domains, y=conservation, ax=ax2, color='crimson', marker='o', linewidth=2.5)
        
        # Formatting
        ax1.set_ylabel('dN/dS Ratio', fontsize=12)
        ax2.set_ylabel('Residue Conservation (%)', fontsize=12)
        ax1.set_xlabel('Protein Domain', fontsize=14)
        ax1.set_title('Evolutionary Selection vs. Functional Conservation', fontsize=16)
        
        # Add neutral evolution line
        plt.axhline(y=1, color='red', linestyle='--', label='Neutral Evolution')
        
        # Add legend
        fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.85))
        
        # Save plot
        plt.tight_layout()
        output_file = output_dir / 'fig1_selection_conservation.png'
        plt.savefig(output_file, dpi=300)
        plt.close()
        
        logging.info(f"Plot saved to: {output_file}")
        
    except Exception as e:
        logging.error(f"Error creating plot: {str(e)}")
        raise

def main():
    try:
        setup_logging()
        logging.info("Starting selection vs conservation plot...")
        plot_selection_conservation()
        logging.info("Plot completed successfully")
        
    except Exception as e:
        logging.error(f"Error in main: {str(e)}")
        raise

if __name__ == "__main__":
    main()

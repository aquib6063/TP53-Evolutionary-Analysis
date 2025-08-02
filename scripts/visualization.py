#!/usr/bin/env python3
"""
TP53 Evolutionary Analysis - Visualization Module
------------------------------------------------
Generates publication-quality figures from analysis results.
"""
import os
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Optional

# Set plotting style
plt.style.use('seaborn')
sns.set_context('paper', font_scale=1.4)

class TP53Visualizer:
    def __init__(self, config_path: str = '../config/analysis_config.yaml'):
        """Initialize the visualizer with configuration."""
        self.config = self._load_config(config_path)
        self.results_dir = Path('../results')
        self.fig_dir = self.results_dir / 'figures'
        self.fig_dir.mkdir(exist_ok=True)
        
        # Set up color palette
        self.palette = sns.color_palette('viridis', n_colors=15)
        
    def _load_config(self, config_path: str) -> Dict:
        """Load configuration from YAML file."""
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    
    def plot_gc_content(self, save: bool = True) -> plt.Figure:
        """Plot GC content comparison across species."""
        # Load data
        df = pd.read_csv(self.results_dir / 'composition_analysis.csv')
        
        # Create figure
        fig, ax = plt.subplots(figsize=(self.config['output']['width'], 
                                      self.config['output']['height']))
        
        # Sort by CTD GC content
        df = df.sort_values('ctd_gc', ascending=False)
        
        # Plot bars
        x = range(len(df))
        width = 0.35
        ax.bar([i - width/2 for i in x], df['cds_gc'], width, label='CDS')
        ax.bar([i + width/2 for i in x], df['ctd_gc'], width, label='CTD')
        
        # Add labels and title
        ax.set_xticks(x)
        ax.set_xticklabels(df['species'], rotation=45, ha='right')
        ax.set_ylabel('GC Content (%)')
        ax.set_title('GC Content Comparison: CDS vs CTD')
        ax.legend()
        
        plt.tight_layout()
        
        if save:
            fig.savefig(self.fig_dir / 'gc_comparison.png', 
                       dpi=self.config['output']['dpi'],
                       bbox_inches='tight')
        
        return fig
    
    def plot_ptm_conservation(self, save: bool = True) -> plt.Figure:
        """Plot PTM site conservation across species."""
        # Load data
        df = pd.read_csv(self.results_dir / 'ptm_analysis.csv')
        
        # Normalize by protein length
        for ptm in ['phosphorylation', 'acetylation', 'ubiquitination']:
            df[f'{ptm}_norm'] = df[ptm] / df['protein_length'] * 100
        
        # Create figure
        fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
        ptm_types = ['phosphorylation', 'acetylation', 'ubiquitination']
        
        for ax, ptm in zip(axes, ptm_types):
            # Sort by normalized count
            df_sorted = df.sort_values(f'{ptm}_norm', ascending=False)
            
            # Plot
            sns.barplot(x='species', y=f'{ptm}_norm', data=df_sorted, 
                       ax=ax, palette=self.palette)
            ax.set_ylabel(f'{ptm.capitalize()} Sites\n(per 100 aa)')
            ax.tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        
        if save:
            fig.savefig(self.fig_dir / 'ptm_conservation.png', 
                       dpi=self.config['output']['dpi'],
                       bbox_inches='tight')
        
        return fig
    
    def plot_phylogeny(self, save: bool = True) -> plt.Figure:
        """Plot phylogenetic tree with bootstrap support."""
        from ete3 import Tree, TreeStyle, NodeStyle, TextFace
        
        # Load tree
        tree_path = self.results_dir / 'phylogeny' / 'tp53.treefile'
        if not tree_path.exists():
            print("Phylogenetic tree not found. Run the phylogeny analysis first.")
            return None
            
        t = Tree(str(tree_path))
        
        # Configure tree style
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.show_branch_support = True
        ts.scale = 100  # Scale branch lengths
        
        # Create figure
        fig = plt.figure(figsize=(self.config['output']['width'], 
                                self.config['output']['height']))
        
        # Render tree
        t.render(str(self.fig_dir / 'phylogeny.png'), tree_style=ts)
        
        return fig

def main():
    """Generate all figures."""
    visualizer = TP53Visualizer()
    
    print("Generating GC content plot...")
    visualizer.plot_gc_content()
    
    print("Generating PTM conservation plot...")
    visualizer.plot_ptm_conservation()
    
    print("Generating phylogenetic tree...")
    visualizer.plot_phylogeny()
    
    print(f"\nAll figures saved to: {visualizer.fig_dir}")

if __name__ == "__main__":
    main()

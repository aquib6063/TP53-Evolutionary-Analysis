import plotly.graph_objects as go
import pandas as pd
import numpy as np
import logging
from pathlib import Path

# Setup logging
def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def create_domain_comparison_plot():
    """Create interactive plot comparing domains across species"""
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
        fig = go.Figure()
        
        # Add GC content bars
        fig.add_trace(go.Bar(
            x=domains,
            y=human_gc,
            name='Human GC%',
            marker_color='blue',
            opacity=0.7
        ))
        
        fig.add_trace(go.Bar(
            x=domains,
            y=pig_gc,
            name='Pig GC%',
            marker_color='orange',
            opacity=0.7
        ))
        
        # Add conservation line
        fig.add_trace(go.Scatter(
            x=domains,
            y=conservation,
            name='Conservation %',
            mode='lines+markers',
            line=dict(color='red', width=3),
            marker=dict(size=10)
        ))
        
        # Update layout
        fig.update_layout(
            title='Domain-wise Comparison of GC Content and Conservation',
            xaxis_title='Protein Domain',
            yaxis_title='Percentage',
            yaxis2=dict(
                title='Conservation %',
                overlaying='y',
                side='right'
            ),
            barmode='group',
            legend_title='Species/Metric',
            height=600,
            width=1000
        )
        
        # Save as HTML
        output_file = output_dir / 'domain_comparison.html'
        fig.write_html(output_file)
        logging.info(f"Interactive plot saved to: {output_file}")
        
    except Exception as e:
        logging.error(f"Error creating domain comparison plot: {str(e)}")
        raise

def create_variant_position_plot():
    """Create interactive plot showing variant positions"""
    try:
        # Create output directory
        output_dir = Path('results/interactive_plots')
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Load variant data
        variants = pd.read_csv('results/variant_analysis/high_impact_variants.csv')
        
        # Create figure
        fig = go.Figure()
        
        # Add scatter plot
        fig.add_trace(go.Scatter(
            x=variants['Position'],
            y=[1] * len(variants),
            mode='markers',
            marker=dict(
                size=10,
                color=variants['IMPACT'].map({
                    'HIGH': 'red',
                    'MODERATE': 'orange',
                    'LOW': 'green'
                }),
                opacity=0.7
            ),
            hovertext=variants['Position'].astype(str) + ': ' + 
                     variants['IMPACT'] + ' - ' + 
                     variants['VARIANT_TYPE'],
            hoverinfo='text'
        ))
        
        # Add domain boundaries
        domain_boundaries = {
            'N-terminal': (1, 42),
            'DNA-binding': (94, 289),
            'Oligomerization': (323, 355),
            'C-terminal': (363, 393)
        }
        
        for name, (start, end) in domain_boundaries.items():
            fig.add_shape(
                type='rect',
                x0=start,
                x1=end,
                y0=0,
                y1=2,
                fillcolor='lightgray',
                opacity=0.3,
                line=dict(width=0)
            )
            
            fig.add_annotation(
                x=(start + end) / 2,
                y=1.5,
                text=name,
                showarrow=False,
                font=dict(size=12)
            )
        
        # Update layout
        fig.update_layout(
            title='Variant Positions Across TP53 Domains',
            xaxis_title='Amino Acid Position',
            yaxis_title='Impact Level',
            yaxis=dict(
                tickvals=[0, 1, 2],
                ticktext=['Low', 'Medium', 'High'],
                range=[-0.5, 2.5]
            ),
            height=600,
            width=1200,
            showlegend=False
        )
        
        # Save as HTML
        output_file = output_dir / 'variant_positions.html'
        fig.write_html(output_file)
        logging.info(f"Interactive plot saved to: {output_file}")
        
    except Exception as e:
        logging.error(f"Error creating variant position plot: {str(e)}")
        raise

def create_dnds_plot():
    """Create interactive plot showing dN/dS ratios"""
    try:
        # Create output directory
        output_dir = Path('results/interactive_plots')
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Domain data
        domains = ['N-terminal', 'DNA-binding', 'Oligomerization', 'C-terminal']
        dnds = [1.2, 0.8, 1.5, 8.49]
        
        # Create figure
        fig = go.Figure()
        
        # Add bars
        fig.add_trace(go.Bar(
            x=domains,
            y=dnds,
            marker_color='skyblue',
            opacity=0.7
        ))
        
        # Add neutral line
        fig.add_shape(
            type='line',
            x0=-0.5,
            x1=3.5,
            y0=1,
            y1=1,
            line=dict(color='red', dash='dash')
        )
        
        # Add annotations
        fig.add_annotation(
            x=3.5,
            y=1,
            text='Neutral Evolution',
            showarrow=False,
            xshift=10,
            font=dict(color='red')
        )
        
        # Update layout
        fig.update_layout(
            title='dN/dS Ratios Across TP53 Domains',
            xaxis_title='Protein Domain',
            yaxis_title='dN/dS Ratio',
            height=600,
            width=1000
        )
        
        # Save as HTML
        output_file = output_dir / 'dnds_ratios.html'
        fig.write_html(output_file)
        logging.info(f"Interactive plot saved to: {output_file}")
        
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

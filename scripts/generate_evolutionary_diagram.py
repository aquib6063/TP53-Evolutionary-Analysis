import matplotlib.pyplot as plt
import networkx as nx
import os

def generate_evolutionary_diagram():
    """Generate evolutionary strategy diagram using matplotlib"""
    # Create graph
    G = nx.DiGraph()
    
    # Add nodes with positions
    pos = {
        'Complex Regulatory Networks': (0, 2),
        'Higher Mutation Risk': (2, 2),
        'Human Cancer': (4, 2),
        'Streamlined Control Mechanisms': (0, 0),
        'Reduced Complexity': (2, 0),
        'Pig Cancer Resistance': (4, 0)
    }
    
    # Add edges
    edges = [
        ('Complex Regulatory Networks', 'Higher Mutation Risk'),
        ('Higher Mutation Risk', 'Human Cancer'),
        ('Streamlined Control Mechanisms', 'Reduced Complexity'),
        ('Reduced Complexity', 'Pig Cancer Resistance')
    ]
    
    G.add_edges_from(edges)
    
    # Create figure
    plt.figure(figsize=(12, 6))
    
    # Draw nodes
    nx.draw_networkx_nodes(G, pos, node_size=0, alpha=0)  # Invisible nodes
    
    # Draw edges
    nx.draw_networkx_edges(G, pos, arrows=True, arrowstyle='->',
                          arrowsize=20, width=2)
    
    # Draw labels
    nx.draw_networkx_labels(G, pos, font_size=12, font_weight='bold')
    
    # Add subgraph labels
    plt.text(-1, 2.5, 'Human Evolution', fontsize=14, fontweight='bold')
    plt.text(-1, 0.5, 'Pig Evolution', fontsize=14, fontweight='bold')
    
    # Add background colors
    plt.fill_between([0, 4], [1.5, 1.5], [2.5, 2.5], color='#f9f9f9', alpha=0.5)
    plt.fill_between([0, 4], [-0.5, -0.5], [0.5, 0.5], color='#e6f0ff', alpha=0.5)
    
    # Remove axes
    plt.axis('off')
    
    # Save figure
    output_file = "results/figures/evolutionary_strategy.png"
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Diagram generated successfully: {output_file}")

if __name__ == "__main__":
    generate_evolutionary_diagram()

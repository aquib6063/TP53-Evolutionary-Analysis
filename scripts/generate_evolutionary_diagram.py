import matplotlib.pyplot as plt
import networkx as nx
import os

def generate_evolutionary_diagram():
    """Generate evolutionary strategy diagram using matplotlib"""
    # Create graph
    G = nx.DiGraph()
    
    # Add nodes with precise positions
    pos = {
        'Complex Regulatory Networks': (-1, 2.5),
        'Higher Mutation Risk': (2, 2.5),
        'Human Cancer': (5, 2.5),
        'Streamlined Control Mechanisms': (-1, -0.5),
        'Reduced Complexity': (2, -0.5),
        'Pig Cancer Resistance': (5, -0.5)
    }
    
    # Add edges with precise positions
    edges = [
        ('Complex Regulatory Networks', 'Higher Mutation Risk'),
        ('Higher Mutation Risk', 'Human Cancer'),
        ('Streamlined Control Mechanisms', 'Reduced Complexity'),
        ('Reduced Complexity', 'Pig Cancer Resistance')
    ]
    
    G.add_edges_from(edges)
    
    # Create figure
    plt.figure(figsize=(12, 6))
    
    # Draw edges with better styling
    for edge in G.edges():
        source, target = edge
        x_values = [pos[source][0], pos[target][0]]
        y_values = [pos[source][1], pos[target][1]]
        plt.plot(x_values, y_values, '->', color='black', linewidth=2, 
                solid_capstyle='round')
        
        # Add arrowhead
        plt.arrow(x_values[1], y_values[1], 0, 0, 
                 head_width=0.2, head_length=0.3, 
                 fc='black', ec='black')
    
    # Draw nodes as rectangles with better alignment
    node_width = 2.5
    node_height = 0.8
    
    for node, (x, y) in pos.items():
        # Draw rectangle
        rect = plt.Rectangle((x - node_width/2, y - node_height/2), 
                           node_width, node_height, 
                           fill=True, edgecolor='black', 
                           facecolor='white', alpha=0.8)
        plt.gca().add_patch(rect)
        
        # Draw text with better alignment
        plt.text(x, y, node, 
                ha='center', va='center', 
                fontsize=12, fontweight='bold', 
                bbox=dict(facecolor='white', alpha=0.0))
    
    # Add subgraph labels with better styling
    plt.text(-2.5, 2.5, 'Human Evolution', 
            fontsize=14, fontweight='bold', 
            bbox=dict(facecolor='#f9f9f9', alpha=0.7))
    
    plt.text(-2.5, -0.5, 'Pig Evolution', 
            fontsize=14, fontweight='bold', 
            bbox=dict(facecolor='#e6f0ff', alpha=0.7))
    
    # Add background colors
    plt.fill_between([-2, 6], [1.5, 1.5], [3.5, 3.5], 
                    color='#f9f9f9', alpha=0.3)
    plt.fill_between([-2, 6], [-1.5, -1.5], [0.5, 0.5], 
                    color='#e6f0ff', alpha=0.3)
    
    # Improve layout
    plt.xlim(-3, 6)
    plt.ylim(-2, 4)
    plt.axis('off')
    plt.tight_layout()
    
    # Save figure with higher quality
    output_file = "results/figures/evolutionary_strategy.png"
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    plt.savefig(output_file, dpi=600, bbox_inches='tight', 
               transparent=True)
    plt.close()
    
    print(f"Diagram generated successfully: {output_file}")

if __name__ == "__main__":
    generate_evolutionary_diagram()

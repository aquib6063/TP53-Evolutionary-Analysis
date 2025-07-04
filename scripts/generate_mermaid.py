import subprocess
import os

def generate_mermaid_diagram():
    """Generate mermaid diagram and save as PNG"""
    # Mermaid diagram code
    mermaid_code = '''
    graph LR
        A[Complex Regulatory Networks] --> B(Higher Mutation Risk) --> C(Human Cancer)
        D[Streamlined Control Mechanisms] --> E(Reduced Complexity) --> F(Pig Cancer Resistance)
        
        style A fill:#f9f,stroke:#333,stroke-width:2px
        style B fill:#f9f,stroke:#333,stroke-width:2px
        style C fill:#f9f,stroke:#333,stroke-width:2px
        style D fill:#bbf,stroke:#333,stroke-width:2px
        style E fill:#bbf,stroke:#333,stroke-width:2px
        style F fill:#bbf,stroke:#333,stroke-width:2px

        subgraph Human Evolution
            A --> B --> C
        end

        subgraph Pig Evolution
            D --> E --> F
        end
    '''
    
    # Save mermaid code to file
    mermaid_file = "mermaid.mmd"
    with open(mermaid_file, 'w') as f:
        f.write(mermaid_code)
    
    # Generate PNG
    output_file = "results/figures/evolutionary_strategy.png"
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Run mermaid-cli
    cmd = f"mermaid {mermaid_file} -o {output_file}"
    subprocess.run(cmd, shell=True, check=True)
    
    # Clean up
    os.remove(mermaid_file)
    
    print(f"Diagram generated successfully: {output_file}")

if __name__ == "__main__":
    generate_mermaid_diagram()

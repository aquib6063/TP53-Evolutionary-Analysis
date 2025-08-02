#!/usr/bin/env python3
import argparse
from Bio import SeqIO
import os

def segment_domains(input_file, output_dir):
    """Segment TP53 domains from aligned sequences."""
    os.makedirs(output_dir, exist_ok=True)
    print(f"Segmenting domains from {input_file}")
    
    # Example domain coordinates (adjust based on your TP53 reference)
    domains = {
        'NTD': (1, 100),     # N-terminal domain
        'DBD': (100, 300),   # DNA-binding domain
        'OD': (300, 350),    # Oligomerization domain
        'CTD': (350, 393)    # C-terminal domain
    }
    
    for record in SeqIO.parse(input_file, "fasta"):
        for domain, (start, end) in domains.items():
            domain_seq = record.seq[start-1:end]  # Convert to 0-based
            output_file = f"{output_dir}/{domain}.fasta"
            with open(output_file, "a") as f:
                f.write(f">{record.id}_{domain}\n{domain_seq}\n")

def main():
    parser = argparse.ArgumentParser(description='Segment TP53 domains')
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    args = parser.parse_args()
    segment_domains(args.input, args.output)

if __name__ == "__main__":
    main()

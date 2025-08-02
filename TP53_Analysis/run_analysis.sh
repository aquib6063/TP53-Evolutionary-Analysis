#!/bin/bash

# Create necessary directories
mkdir -p alignments/domains results/ptm_analysis logs

echo "=== TP53 Analysis Pipeline ==="
echo "1. Running multiple sequence alignment..."
mafft --auto --adjustdirection combined_sequences.fasta > alignments/tp53_aligned.fasta

echo "2. Building phylogenetic tree..."
iqtree -s alignments/tp53_aligned.fasta -m MFP -bb 1000 -nt AUTO -pre results/tp53_tree

echo "3. Segmenting domains..."
python3 scripts/segment_domains.py alignments/tp53_aligned.fasta alignments/domains

echo "4. Running base composition analysis..."
python3 scripts/analyze_composition.py alignments/tp53_aligned.fasta results/base_composition.csv

echo -e "\n=== Analysis Complete ==="
echo "Results are saved in the 'results' directory"

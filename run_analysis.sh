#!/bin/bash
# TP53 Evolutionary Analysis Pipeline

# Set up directories
BASE_DIR=$(pwd)
RAW_DATA="$BASE_DIR/data/raw/sequences"
PROCESSED_DIR="$BASE_DIR/data/processed"
RESULTS_DIR="$BASE_DIR/results"
SCRIPTS_DIR="$BASE_DIR/scripts"

echo "Starting TP53 Analysis"
echo "====================="
echo "Base directory: $BASE_DIR"
echo "Raw data: $RAW_DATA"
echo "Processed data: $PROCESSED_DIR"
echo "Results: $RESULTS_DIR"

# Create output directories
mkdir -p "$PROCESSED_DIR/alignments/domains" \
         "$RESULTS_DIR/trees" \
         "$RESULTS_DIR/alignments" \
         "$RESULTS_DIR/ptm_analysis"

# 1. Run alignment
echo -e "\n1. Running multiple sequence alignment..."
mafft --auto "$RAW_DATA/combined.fasta" > "$PROCESSED_DIR/alignments/tp53_aligned.fasta"

# 2. Build tree
echo -e "\n2. Building phylogenetic tree..."
iqtree -s "$PROCESSED_DIR/alignments/tp53_aligned.fasta" \
       -m MFP \
       -bb 1000 \
       -nt AUTO \
       -pre "$RESULTS_DIR/trees/tp53"

echo -e "\nAnalysis complete! Results are in the 'results' directory."

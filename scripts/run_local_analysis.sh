#!/bin/bash

# Check if Java is installed
if ! command -v java &> /dev/null; then
    echo "Java is not installed. Please install Java first."
    echo "You can install Java using Homebrew:"
    echo "brew install openjdk"
    exit 1
fi

# Check if SnpEff is installed
if [ ! -f "snpEff/snpEff.jar" ]; then
    echo "Downloading SnpEff..."
    wget -O snpEff.zip https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
    unzip snpEff.zip
fi

# Convert variants to VCF format
python scripts/create_vcf.py results/variation_analysis/pig_specific_variations.csv results/variation_analysis/pig_variants.vcf

# Run SnpEff annotation
java -jar snpEff/snpEff.jar Sscrofa11.1 results/variation_analysis/pig_variants.vcf > results/variation_analysis/annotated_variants.vcf

# Analyze the variants
python scripts/analyze_variants.py

#!/bin/bash

# Build the Docker image
docker build -t tp53-variant-analysis .

# Run the analysis in Docker
docker run -v $(pwd):/data tp53-variant-analysis bash -c "
    # Copy data to container
    cp /data/results/variation_analysis/pig_variants.vcf /opt/analysis_scripts/
    
    # Run SnpEff annotation
    java -jar /opt/snpEff/snpEff.jar Sscrofa11.1 /opt/analysis_scripts/pig_variants.vcf > /opt/analysis_scripts/annotated_variants.vcf
    
    # Run VEP annotation
    perl /opt/vep/variant_effect_predictor/variant_effect_predictor.pl \
        -i /opt/analysis_scripts/annotated_variants.vcf \
        -o /opt/analysis_scripts/vep_annotated.txt \
        --species pig \
        --cache
    
    # Run variant analysis
    python3 /opt/analysis_scripts/analyze_variants.py
    
    # Copy results back to host
    cp /opt/analysis_scripts/results/* /data/results/
"

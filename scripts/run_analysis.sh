#!/bin/bash

# Run variant analysis
python scripts/analyze_pig_variants.py

# Check if analysis was successful
if [ $? -eq 0 ]; then
    echo "Analysis completed successfully!"
    echo "Results saved in results/variant_analysis/"
else
    echo "Error in analysis. Please check the logs."
fi

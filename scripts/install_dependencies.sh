#!/bin/bash

# Install Python packages
pip install biopython pandas requests beautifulsoup4

# Check if installation was successful
if [ $? -eq 0 ]; then
    echo "Dependencies installed successfully!"
else
    echo "Error installing dependencies"
    exit 1
fi

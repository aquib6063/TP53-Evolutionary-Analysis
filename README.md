# TP53 Evolutionary Analysis

This repository contains comprehensive analysis code and results for studying TP53 evolutionary patterns in pigs and other species, focusing on oncogenic enrichment and functional conservation.

## Project Overview

This project investigates the evolutionary patterns of TP53 across multiple species, with a focus on identifying oncogenic enrichment patterns and functional conservation in pig TP53. The analysis includes:
│   ├── 3_analyze_variations.py
│   ├── 4_gc_content_analysis.py
│   └── 5_generate_figures.py
├── results/
│   ├── variation_analysis/
│   │   ├── pig_specific_variants.csv
│   │   └── domain_variation_summary.csv
│   ├── gc_content/
│   │   └── gc_by_domain.csv
│   └── figures/
│       ├── alignment_conservation.png
│       ├── gc_content_changes.png
│       └── domain_variation_plot.png
├── docs/
│   ├── METHODS.md
│   └── INTERPRETATION.md
├── environment.yml
├── requirements.txt
└── README.md
```

## Setup

1. Create a conda environment:
```bash
conda env create -f environment.yml
conda activate tp53_analysis
```

2. Run analysis pipeline:
```bash
# Run scripts in order
python scripts/1_download_sequences.py
python scripts/2_run_alignment.py
python scripts/3_analyze_variations.py
python scripts/4_gc_content_analysis.py
python scripts/5_generate_figures.py
python scripts/6_combine_fasta.py
```

## Data

The project uses NCBI datasets for TP53 sequences across five species:
- Human
- Pig
- Rhesus Monkey
- Zebrafish
- Rat

## Results

The analysis generates:
1. Sequence alignment
2. Variation analysis across species
3. GC content analysis
4. Visualization plots

## Documentation

Detailed documentation is available in the `docs/` directory:
- `METHODS.md`: Analysis methodology
- `INTERPRETATION.md`: Results interpretation

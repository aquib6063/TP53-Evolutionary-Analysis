w# TP53 Evolutionary Analysis Pipeline

This repository contains a comprehensive pipeline for analyzing TP53 evolution across multiple species, with a focus on the C-terminal domain (CTD) and its role in cancer resistance.

## Table of Contents
- [Installation](#installation)
- [Setup](#setup)
- [Running the Analysis](#running-the-analysis)
- [Output Files](#output-files)
- [Analysis Components](#analysis-components)
- [Troubleshooting](#troubleshooting)

## Installation

### Prerequisites
- Python 3.8+
- Conda (recommended)
- Required bioinformatics tools:
  - MAFFT (for sequence alignment)
  - IQ-TREE (for phylogenetic analysis)
  - DSSP (for secondary structure prediction)
  - IUPred (for disorder prediction)

### Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/TP53-Evolutionary-Analysis.git
   cd TP53-Evolutionary-Analysis
   ```

2. Create and activate the conda environment:
   ```bash
   conda env create -f environment.yml
   conda activate tp53-env
   ```

3. Install Python dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Running the Analysis

### 1. Prepare Input Data
Place your FASTA files in the `data/raw/` directory with the naming convention: `{species}.fasta`

### 2. Run the Complete Pipeline
```bash
python scripts/run_analysis.py
```

### 3. Run Individual Components

#### Base Composition Analysis
```bash
python scripts/run_analysis.py --analysis composition
```

#### Phylogenetic Analysis
```bash
python scripts/run_analysis.py --analysis phylogeny
```

#### PTM Prediction
```bash
python scripts/run_analysis.py --analysis ptm
```

## Output Files

All results are saved in the `results/` directory:

- `composition/`: Base composition analysis
- `phylogeny/`: Phylogenetic trees and alignments
- `ptm/`: PTM prediction results
- `structure/`: Structural analysis results
- `tables/`: Summary tables in CSV format
- `figures/`: Publication-quality figures

## Analysis Components

### 1. Base Composition Analysis
- GC content (CDS and CTD)
- Nucleotide frequencies
- GC skew analysis
- CpG island prediction

### 2. Evolutionary Analysis
- Multiple sequence alignment
- dN/dS calculation
- Positive selection tests
- Phylogenetic tree construction

### 3. PTM Prediction
- Phosphorylation sites
- Acetylation sites
- Ubiquitination sites
- Conservation analysis

### 4. Structural Analysis
- Secondary structure prediction
- Disorder prediction
- Surface accessibility
- Structural alignment

## Troubleshooting

### Common Issues
1. **Missing Dependencies**:
   ```bash
   conda install -c bioconda mafft iqtree
   ```

2. **Memory Issues**:
   - Reduce the number of threads in `config/analysis_config.yaml`
   - Increase system swap space if needed

3. **File Permissions**:
   ```bash
   chmod +x scripts/*.py
   ```

## Citation
If you use this pipeline in your research, please cite:
```
[Citation will be added upon publication]
```

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

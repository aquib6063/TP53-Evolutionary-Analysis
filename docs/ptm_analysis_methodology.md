# TP53 PTM Analysis Methodology

## Overview
This document outlines the methodology used for analyzing post-translational modification (PTM) disruptions between human and pig TP53 sequences. The analysis focuses on identifying and characterizing PTM site differences that could impact protein function and regulation.

## Data Sources

### Sequence Data
- Human TP53: UniProt ID P04637
- Pig TP53: UniProt ID A0A287BKT5

Sequences are retrieved from UniProt using the ExPASy API.

## PTM Prediction Methodology

### PTM Residue Types
The analysis considers the following PTM types based on residue-specific modifications:

| Residue | Possible PTM Types |
|---------|-------------------|
| S       | phosphorylation   |
| T       | phosphorylation   |
| Y       | phosphorylation   |
| K       | acetylation, methylation, ubiquitination |
| R       | methylation       |
| E       | ubiquitination    |
| N       | ubiquitination    |
| C       | ubiquitination    |

### Prediction Scoring
The PTM prediction uses a combination of scoring methods:

1. **Local Context Analysis**
   - Analyzes 5-residue windows around each PTM candidate
   - Calculates conservation scores based on residue context
   - Considers sequence similarity and residue environment

2. **Thresholds**
   - Phosphorylation: 0.5
   - Acetylation: 0.4
   - Methylation: 0.3
   - Ubiquitination: 0.2

## Domain Analysis
The analysis is performed on multiple domains to capture regional differences:

1. Domain 1: 360-393
2. Domain 2: 350-400
3. Domain 3: 325-393
4. Domain 4: 300-450
5. Domain 5: 330-380
6. Domain 6: 340-410

## Analysis Pipeline

### 1. Sequence Retrieval
- Fetch sequences from UniProt
- Validate sequence lengths
- Perform pairwise alignment

### 2. PTM Prediction
- Apply local context analysis
- Calculate conservation scores
- Generate PTM site predictions

### 3. Disruption Analysis
- Compare PTM sites between species
- Identify residue changes
- Classify disruption types:
  - Loss: PTM site present in human but missing in pig
  - Gain: PTM site present in pig but missing in human
  - Change: PTM type change between species

### 4. Statistical Analysis
- Calculate sequence similarity scores
- Analyze PTM type distributions
- Identify common disruption patterns

## Visualization

### PTM Type Distribution
- Bar plots showing PTM type frequencies
- Separate plots for human and pig sequences
- Comparison of PTM type distributions

### PTM Changes
- Heatmap showing PTM type transitions
- Visualization of disruption patterns
- Domain-wise analysis plots

## Results Interpretation

### Key Metrics
- Total disruptions per domain
- Sequence similarity scores
- Common PTM changes
- Residue-specific modifications

### Impact Analysis
- Functional implications of PTM changes
- Potential impact on protein regulation
- Conservation analysis across domains

## Output Files

1. **CSV Files**
   - `ptm_disruptions_combined.csv`: Combined results
   - `ptm_disruptions_domain_*.csv`: Domain-specific results
   - `domain_analysis.csv`: Domain-wise statistics

2. **Visualization Files**
   - `ptm_type_distribution.png`: PTM type comparison
   - `ptm_type_changes.png`: PTM change heatmap
   - `domain_analysis.png`: Domain-wise analysis

3. **Log Files**
   - `ptm_analysis.log`: Detailed analysis logs
   - `ptm_type_analysis.log`: PTM type analysis logs

## Future Work
- Integration with additional PTM databases
- Machine learning-based prediction improvement
- Functional impact analysis
- Cross-species comparison expansion

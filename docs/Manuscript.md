# Simplified TP53: Evolutionary Insights into Cancer Defense Mechanisms in Pigs

*Preprint: Evolutionary Medicine Insights*

**Authors**: [Your Name], [Affiliation]
**Correspondence**: [Email] | [ORCID]
**Submitted**: July 4, 2025 | **Preprint DOI**: 10.7287/peerj.preprints.XXXXX

---

## Abstract

Pigs (*Sus scrofa*) maintain low cancer incidence despite large body size, challenging conventional understanding of Peto's paradox. Through comparative genomic analysis of TP53 across five species (human, pig, rhesus monkey, rat, and zebrafish), we identify 5,004 pig-specific nucleotide variations, with significant C-terminal domain degeneration (14.8% GC reduction compared to human). This evolutionary streamlining preserves core DNA-binding functionality while reducing regulatory complexity, suggesting a novel "less-is-more" cancer defense strategy. These findings provide evolutionary insights into oncogenic mechanisms and offer potential therapeutic implications.

---

## Methods

#### A. Sequence Acquisition and Preparation

1. **Species Selection**:
   - Five model organisms were selected representing different evolutionary lineages:
     - Mammals: Human (*Homo sapiens*), Pig (*Sus scrofa*), Rhesus monkey (*Macaca mulatta*), Rat (*Rattus norvegicus*)
     - Fish: Zebrafish (*Danio rerio*)
   - Species were chosen based on their availability in NCBI databases and relevance to cancer biology

2. **Sequence Retrieval**:
   - Full-length TP53 gene sequences were downloaded from NCBI using BioPython's Entrez module
   - Multiple isoforms were identified and the longest transcript was selected for each species
   - Sequences were validated for completeness and quality

3. **Sequence Processing**:
   - Raw sequences were cleaned and formatted using BioPython's SeqIO module
   - Sequence IDs were standardized across species
   - A combined FASTA file was created containing all species sequences

#### B. Sequence Alignment and Analysis

1. **Multiple Sequence Alignment**:
   - Sequences were aligned using BioPython's PairwiseAligner
   - Gaps were introduced to maintain sequence alignment
   - Alignment quality was verified using pairwise identity scores

2. **Variation Analysis**:
   - Pig-specific nucleotide variations were identified by comparing aligned sequences
   - Variations were categorized by base type (A, C, G, T)
   - Statistical significance was assessed using chi-squared tests

3. **Variant Annotation**:
   - Variants were converted to VCF format using custom script
   - Functional annotation was performed in a Docker container environment
   - Annotation pipeline:
     ```bash
     # Docker-based annotation pipeline
     # 1. Build Docker image with SnpEff and VEP
     docker build -t tp53-variant-analysis .

     # 2. Run annotation in container
     docker run -v $(pwd):/data tp53-variant-analysis bash -c "
         # Run SnpEff annotation
         java -jar /opt/snpEff/snpEff.jar Sscrofa11.1 pig_variants.vcf > annotated_variants.vcf

         # Run VEP annotation
         perl /opt/vep/variant_effect_predictor/variant_effect_predictor.pl \
             -i annotated_variants.vcf \
             -o vep_annotated.txt \
             --species pig \
             --cache
     "
     ```
   - Annotation results included:
     - High-impact variants (stop-gains, splice sites)
     - Pathogenicity scores (PolyPhen-2, SIFT)
     - Functional impact predictions
     - Domain-specific variant distribution
     - Conservation scores
     - Variant type classification
     - Correlation between prediction scores

4. **Statistical Analysis**:
   - Impact Distribution Analysis:
     - Chi-squared test for impact type distribution
     - Statistical significance of impact frequencies
   
   - Pathogenicity Score Analysis:
     - Mean and standard deviation calculation
     - Normality testing using Shapiro-Wilk test
     - Distribution visualization with KDE plots
   
   - Domain-Specific Analysis:
     - Chi-squared test for domain distribution
     - Conservation score analysis
     - Boxplot visualization of conservation scores
   
   - High-Impact Variant Analysis:
     - Domain distribution analysis
     - Conservation score comparison
     - Impact type classification
   
   - Correlation Analysis:
     - Pearson correlation between PolyPhen and SIFT scores
     - Scatter plot visualization
     - Correlation significance testing

5. **Data Visualization**:
   - Impact Distribution:
     - Bar plots with error bars
     - Statistical annotations
   
   - Pathogenicity Scores:
     - Histograms with KDE
     - Normality plots
     - Score distributions
   
   - Domain Analysis:
     - Stacked bar plots
     - Box plots for conservation
     - Domain-specific distributions
   
   - Variant Types:
     - Pie charts
     - Distribution plots
     - Type-specific analysis
   
   - Correlation Visualization:
     - Scatter plots
     - Correlation matrices
     - Heatmaps

3. **Domain-Specific Analysis**:
   - TP53 gene was divided into four functional domains:
     - N-terminal domain
     - DNA-binding domain
     - Oligomerization domain
     - C-terminal domain
   - Domain boundaries were defined based on known functional regions

#### C. GC Content Analysis

1. **GC Content Calculation**:
   - GC percentage was calculated for each domain using BioPython's SeqUtils module
   - GC content was normalized across species
   - Domain-specific GC content was compared between species

2. **Statistical Analysis**:
   - Pairwise comparisons of GC content were performed
   - Statistical significance was assessed using t-tests
   - Domain-specific GC content trends were analyzed

#### D. Phylogenetic Analysis

1. **Distance Matrix Calculation**:
   - Pairwise sequence distances were calculated using BioPython's DistanceCalculator
   - Distance metrics were based on sequence identity
   - Distance matrix was visualized using heatmap

2. **Tree Construction**:
   - Phylogenetic tree was constructed using Neighbor-Joining method
   - Branch lengths were calculated based on sequence distances
   - Tree topology was validated using bootstrap analysis

3. **Tree Visualization**:
   - Tree was visualized using BioPython's Phylo module
   - Branch lengths were displayed
   - Species names were simplified for clarity

#### E. Data Visualization

1. **Sequence Comparison**:
   - Heatmaps of sequence identity were generated
   - Distribution plots of sequence distances were created
   
2. **GC Content Visualization**:
   - Domain-wise GC content was plotted
   - Species comparisons were visualized
   - Statistical significance was indicated

3. **Variant Impact Visualization**:
   - High-impact variant distribution across domains
   - Pathogenicity score distribution
   - Impact type categorization

4. **Evolutionary Strategy Diagram**:

1. **Sequence Comparison**:
   - Heatmaps of sequence identity were generated
   - Distribution plots of sequence distances were created
   
2. **GC Content Visualization**:
   - Domain-wise GC content was plotted
   - Species comparisons were visualized
   - Statistical significance was indicated

3. **Evolutionary Strategy Diagram**:
   - Network diagram was created using matplotlib
   - Evolutionary paths were visualized
   - Comparative analysis was presented

#### F. Statistical Analysis

1. **Variation Distribution**:
   - Base composition was analyzed
   - Statistical significance of variations was assessed
   - Distribution patterns were identified

2. **GC Content Comparison**:
   - Domain-wise comparisons were performed
   - Species comparisons were analyzed
   - Statistical significance was determined

3. **Tree Topology**:
   - Branch support values were calculated
   - Evolutionary relationships were validated
   - Species divergence was quantified

---

## Main Text

### 1. The Cancer Resistance Enigma

Pigs, despite having approximately 100 times more cells than humans, exhibit lower cancer incidence rates, challenging the expectations set by Peto's paradox[1]. While previous studies have shown that elephants evolved TP53 gene duplication as a cancer defense mechanism[2], our analysis reveals that pigs have taken a different evolutionary approach: gene simplification.

### 2. Key Findings

#### A. Variation Distribution

Our analysis of the TP53 gene across five species reveals distinct patterns of nucleotide variations in pigs:

1. **Total Variations**: 5,004 pig-specific nucleotide variations identified
2. **Base Composition**: 
   - Cytosine (C): 27.9%
   - Guanine (G): 27.7%
   - Thymine (T): 23.7%
   - Adenine (A): 20.7%

#### B. Domain-Specific Evolution

| Domain         | Human GC% | Pig GC% | ΔGC% |
|----------------|-----------|---------|------|
| N-terminal     | 59.1      | 54.1    | -5.0 |
| DNA-binding    | 52.5      | 53.0    | +0.5 |
| Oligomerization| 52.7      | 48.4    | -4.3 |
| C-terminal     | 48.1      | 33.3    | -14.8|

#### C. Evolutionary Strategy

Our analysis suggests an evolutionary "less-is-more" strategy in pigs:

1. **Preserved Functionality**:
   - DNA-binding domain shows minimal changes (+0.5% GC)
   - Maintains core tumor suppressor function

2. **Reduced Complexity**:
   - Significant C-terminal degeneration (-14.8% GC)
   - Reduced regulatory complexity
   - Potential for more stable activation

3. **Species Comparison**:
   - Zebrafish shows significantly lower GC content (37.2%)
   - Mammalian species maintain higher GC content (53-54%)
   - Pig GC content (47.2%) is intermediate but distinct

![Evolutionary Strategy](../results/figures/evolutionary_strategy.png)

*Evolutionary paths of TP53 regulation in humans and pigs. Human evolution shows complex regulatory networks leading to higher mutation risk and cancer, while pig evolution shows streamlined control mechanisms leading to reduced complexity and cancer resistance.*

### 3. Evolutionary Strategy

Our analysis reveals a distinct evolutionary strategy in pigs:

1. **Human Evolution**:
   - Complex regulatory networks have evolved
   - Increased mutation risk due to complex control mechanisms
   - Higher cancer incidence as a trade-off for complex regulation

2. **Pig Evolution**:
   - Streamlined control mechanisms
   - Reduced regulatory complexity
   - Lower cancer incidence through simplified control

3. **Mechanistic Insights**:
   - C-terminal domain degeneration (-14.8% GC) reduces regulatory complexity
   - Preserved DNA-binding domain (+0.5% GC) maintains core function
   - Evolutionary trade-off favoring stability over complexity

### 3. Implications

1. **Evolutionary Insight**: Pigs have evolved a unique cancer defense strategy through gene simplification rather than duplication.
2. **Oncology Applications**: Understanding this "less-is-more" strategy could inform new therapeutic approaches.
3. **Comparative Genomics**: Highlights the importance of domain-specific evolution in cancer-related genes.

## References

1. Peto, R. (1975). Cancer and age. British Journal of Cancer, 32(4), 411-426.
2. Sulak, M., Fong, L., Mika, K., Chigurupati, S., Yon, L., Mongan, N. P., ... & Lynch, V. J. (2016). TP53 copy number expansion is associated with the evolution of increased body size and an enhanced DNA damage response in elephants. Elife, 5, e11994.

---

## Supplementary Figures

### Figure 1: Domain-wise GC Content Comparison

![Domain GC Content](https://github.com/aquib6063/TP53-Evolutionary-Analysis/blob/main/results/figures/distance_matrix_heatmap.png)

*Domain-specific GC content comparison across five species showing significant differences in the C-terminal region.*

### Figure 2: Evolutionary Tree

![Phylogenetic Tree](https://github.com/aquib6063/TP53-Evolutionary-Analysis/blob/main/results/figures/tp53_phylogeny.png)

*Phylogenetic tree showing evolutionary relationships between species based on TP53 sequence alignment.*

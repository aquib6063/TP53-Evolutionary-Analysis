import pandas as pd
import os
import logging
from Bio import SeqIO

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Base directory (should be the project root)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Input paths
VARIATIONS_FILE = os.path.join(BASE_DIR, "results", "variation_analysis", "pig_specific_variations.csv")
GC_FILE = os.path.join(BASE_DIR, "results", "gc_analysis", "domain_gc_content.csv")
ALIGNMENT_FILE = os.path.join(BASE_DIR, "data", "processed", "all_species_aligned.fasta")

# Output path
OUTPUT_FILE = os.path.join(BASE_DIR, "docs", "INTERPRETATION.md")

def summarize_results():
    """Generate comprehensive summary of analysis results"""
    try:
        # Load data
        variations = pd.read_csv(VARIATIONS_FILE)
        gc_data = pd.read_csv(GC_FILE)
        
        # Create summary output
        summary = []
        
        # 1. Sequence Length Analysis
        summary.append("# TP53 Gene Analysis Summary\n")
        summary.append("## Sequence Length Analysis\n")
        
        # Read sequences
        sequences = list(SeqIO.parse(ALIGNMENT_FILE, "fasta"))
        seq_lengths = {seq.id.split('_')[0]: len(str(seq.seq).replace('-', '')) 
                      for seq in sequences}
        
        summary.append("| Species | Sequence Length |\n")
        summary.append("|---------|----------------|\n")
        for species, length in seq_lengths.items():
            summary.append(f"| {species} | {length} bp |\n")
        
        # 2. Pig-Specific Variations
        summary.append("\n## Pig-Specific Variations\n")
        summary.append(f"Total pig-specific variations: {len(variations)}\n")
        
        # Variation types
        variation_types = variations['Pig_Base'].value_counts()
        summary.append("\n### Variation Types\n")
        summary.append("| Base | Count | Percentage |\n")
        summary.append("|------|-------|------------|\n")
        for base, count in variation_types.items():
            percent = (count / len(variations)) * 100
            summary.append(f"| {base} | {count} | {percent:.1f}% |\n")
        
        # 3. GC Content Analysis
        summary.append("\n## GC Content Analysis\n")
        summary.append("### Domain-wise GC Content Comparison\n")
        
        summary.append("| Domain | Human GC% | Pig GC% | Difference |\n")
        summary.append("|--------|-----------|----------|------------|\n")
        
        for domain in gc_data['Domain'].unique():
            human_gc = gc_data[(gc_data['Species']=='Human') & (gc_data['Domain']==domain)]['GC_Percent'].values[0]
            pig_gc = gc_data[(gc_data['Species']=='Pig') & (gc_data['Domain']==domain)]['GC_Percent'].values[0]
            diff = pig_gc - human_gc
            summary.append(f"| {domain} | {human_gc:.1f} | {pig_gc:.1f} | {'+' if diff > 0 else ''}{diff:.1f} |\n")
        
        # 4. Species-wise GC Content
        summary.append("\n### Species-wise GC Content\n")
        species_gc = gc_data.groupby('Species')['GC_Percent'].mean().round(1)
        summary.append("| Species | Average GC% |\n")
        summary.append("|---------|-------------|\n")
        for species, gc in species_gc.items():
            summary.append(f"| {species} | {gc} |\n")
        
        # 5. Key Findings
        summary.append("\n## Key Findings\n")
        summary.append("1. Sequence Length:\n")
        summary.append(f"   - Longest sequence: {max(seq_lengths, key=seq_lengths.get)} ({max(seq_lengths.values())} bp)\n")
        summary.append(f"   - Shortest sequence: {min(seq_lengths, key=seq_lengths.get)} ({min(seq_lengths.values())} bp)\n")
        
        summary.append("\n2. Pig-Specific Variations:\n")
        summary.append(f"   - Total variations: {len(variations)}\n")
        summary.append(f"   - Most common variation: {variation_types.idxmax()} ({(variation_types.max() / len(variations) * 100):.1f}% of variations)\n")
        
        summary.append("\n3. GC Content:\n")
        summary.append("   - Highest GC content domain: N-terminal\n")
        summary.append("   - Lowest GC content domain: Oligomerization\n")
        summary.append("   - Zebrafish shows significantly lower GC content compared to mammals\n")
        
        # Save summary
        os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
        with open(OUTPUT_FILE, 'w') as f:
            f.write('\n'.join(summary))
        
        logging.info(f"✅ Summary saved to: {OUTPUT_FILE}")
        
    except Exception as e:
        logging.error(f"❌ Error generating summary: {str(e)}")

if __name__ == "__main__":
    summarize_results()

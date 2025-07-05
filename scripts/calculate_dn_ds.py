import logging
from Bio.Seq import Seq
from Bio.SeqUtils import CodonUsage
from Bio.Data import CodonTable
from pathlib import Path
import pandas as pd

# Setup logging
def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def count_synonymous_sites(codon_alignment):
    """Count synonymous sites in codon alignment"""
    standard_table = CodonTable.unambiguous_dna_by_id[1]
    syn_sites = 0
    
    # Count sites where synonymous mutations are possible
    for codon in codon_alignment:
        amino_acid = standard_table.forward_table.get(codon, '*')
        if amino_acid != '*':
            # Get all possible codons for this amino acid
            possible_codons = [c for c, aa in standard_table.forward_table.items() 
                             if aa == amino_acid]
            
            # Count sites where a single nucleotide change results in same AA
            for i in range(3):
                # Generate all possible single nucleotide changes
                for base in 'ACGT':
                    if base != codon[i]:
                        new_codon = codon[:i] + base + codon[i+1:]
                        if new_codon in possible_codons:
                            syn_sites += 1
                            break
    
    return syn_sites / 3  # Divide by 3 since each codon has 3 sites

def count_nonsynonymous_sites(codon_alignment):
    """Count nonsynonymous sites in codon alignment"""
    standard_table = CodonTable.unambiguous_dna_by_id[1]
    non_syn_sites = 0
    
    for codon in codon_alignment:
        amino_acid = standard_table.forward_table.get(codon, '*')
        if amino_acid != '*':
            # Get all possible codons for this amino acid
            possible_codons = [c for c, aa in standard_table.forward_table.items() 
                             if aa == amino_acid]
            
            # Count sites where a single nucleotide change results in different AA
            for i in range(3):
                # Generate all possible single nucleotide changes
                for base in 'ACGT':
                    if base != codon[i]:
                        new_codon = codon[:i] + base + codon[i+1:]
                        if new_codon not in possible_codons:
                            non_syn_sites += 1
                            break
    
    return non_syn_sites / 3

def count_substitutions(ref_seq, alt_seq):
    """Count substitutions between two sequences"""
    syn_substitutions = 0
    non_syn_substitutions = 0
    
    standard_table = CodonTable.unambiguous_dna_by_id[1]
    
    # Process each codon position
    for i in range(0, len(ref_seq), 3):
        ref_codon = ref_seq[i:i+3]
        alt_codon = alt_seq[i:i+3]
        
        if ref_codon == alt_codon:
            continue
            
        # Count substitutions
        if ref_codon in standard_table.forward_table and \
           alt_codon in standard_table.forward_table:
            
            ref_aa = standard_table.forward_table[ref_codon]
            alt_aa = standard_table.forward_table[alt_codon]
            
            if ref_aa == alt_aa:
                syn_substitutions += 1
            else:
                non_syn_substitutions += 1
    
    return syn_substitutions, non_syn_substitutions

def calculate_dnds(ref_seq, alt_seq):
    """Calculate dN/dS ratio between two sequences"""
    try:
        # Convert to uppercase and remove gaps
        ref_seq = str(ref_seq).upper().replace('-', '')
        alt_seq = str(alt_seq).upper().replace('-', '')
        
        # Ensure sequences are same length and divisible by 3
        if len(ref_seq) != len(alt_seq) or len(ref_seq) % 3 != 0:
            raise ValueError("Sequences must be same length and codon-aligned")
            
        # Count sites
        syn_sites = count_synonymous_sites(ref_seq)
        non_syn_sites = count_nonsynonymous_sites(ref_seq)
        
        # Count substitutions
        syn_substitutions, non_syn_substitutions = count_substitutions(ref_seq, alt_seq)
        
        # Calculate dN and dS
        dN = non_syn_substitutions / non_syn_sites
        dS = syn_substitutions / syn_sites
        
        return dN/dS if dS > 0 else float('inf')
        
    except Exception as e:
        logging.error(f"Error calculating dN/dS: {str(e)}")
        raise

def process_alignment(alignment_file):
    """Process multiple sequence alignment and calculate dN/dS for each pair"""
    try:
        # Create output directory
        output_dir = Path('results/dnds_analysis')
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Read alignment
        alignment = pd.read_csv(alignment_file)
        
        # Calculate dN/dS for each pair
        results = []
        for i in range(len(alignment) - 1):
            for j in range(i + 1, len(alignment)):
                ref_seq = alignment.iloc[i]['sequence']
                alt_seq = alignment.iloc[j]['sequence']
                
                try:
                    dnds = calculate_dnds(ref_seq, alt_seq)
                    results.append({
                        'species1': alignment.iloc[i]['species'],
                        'species2': alignment.iloc[j]['species'],
                        'dN/dS': dnds
                    })
                except Exception as e:
                    logging.error(f"Error processing pair {i}-{j}: {str(e)}")
                    continue
        
        # Save results
        results_df = pd.DataFrame(results)
        output_file = output_dir / 'dnds_results.csv'
        results_df.to_csv(output_file, index=False)
        logging.info(f"dN/dS results saved to: {output_file}")
        
        return results_df
        
    except Exception as e:
        logging.error(f"Error processing alignment: {str(e)}")
        raise

def main():
    try:
        setup_logging()
        logging.info("Starting dN/dS calculation...")
        
        # Process alignment
        alignment_file = 'results/alignment_and_tree/aligned_sequences.csv'
        results = process_alignment(alignment_file)
        
        # Calculate domain-specific dN/dS
        domain_boundaries = {
            'N-terminal': (0, 42*3),  # Convert to nucleotide positions
            'DNA-binding': (94*3, 289*3),
            'Oligomerization': (323*3, 355*3),
            'C-terminal': (363*3, 393*3)
        }
        
        for domain, (start, end) in domain_boundaries.items():
            domain_results = []
            for _, row in results.iterrows():
                ref_seq = row['species1_sequence'][start:end]
                alt_seq = row['species2_sequence'][start:end]
                try:
                    dnds = calculate_dnds(ref_seq, alt_seq)
                    domain_results.append({
                        'species1': row['species1'],
                        'species2': row['species2'],
                        'domain': domain,
                        'dN/dS': dnds
                    })
                except Exception as e:
                    logging.error(f"Error processing {domain} domain: {str(e)}")
                    continue
            
            domain_df = pd.DataFrame(domain_results)
            output_file = f'results/dnds_analysis/{domain}_dnds.csv'
            domain_df.to_csv(output_file, index=False)
            logging.info(f"Domain-specific dN/dS saved to: {output_file}")
        
        logging.info("dN/dS calculation completed successfully")
        
    except Exception as e:
        logging.error(f"Error in main: {str(e)}")
        raise

if __name__ == "__main__":
    main()

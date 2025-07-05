import os
import logging
from pathlib import Path
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
import json
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from Bio.Phylo.PAML import codeml
import datetime

def setup_logging():
    """Setup logging configuration"""
    logging.basicConfig(
        filename='selection_analysis.log',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def calculate_dn_ds(alignment_file, tree_file):
    """Calculate dN/dS ratio using BioPython"""
    try:
        # Read alignment
        alignment = AlignIO.read(alignment_file, 'phylip-relaxed')
        
        # Calculate distances
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(alignment)
        
        # Calculate dN/dS
        # Convert alignment to codons
        codon_alignment = MultipleSeqAlignment([
            SeqRecord(seq.seq, id=seq.id)
            for seq in alignment
        ])
        
        # Calculate number of non-synonymous and synonymous sites
        non_synonymous = 0
        synonymous = 0
        
        for i in range(0, len(codon_alignment[0].seq), 3):
            codons = [str(seq.seq[i:i+3]) for seq in codon_alignment]
            if '-' in codons[0] or '-' in codons[1]:
                continue
            
            # Skip if any codon is invalid (e.g., stop codon)
            if any(codon not in CodonTable.unambiguous_dna_by_name["Standard"].forward_table for codon in codons):
                continue
            
            if codons[0] != codons[1]:
                if CodonTable.unambiguous_dna_by_name["Standard"].forward_table[codons[0]] != CodonTable.unambiguous_dna_by_name["Standard"].forward_table[codons[1]]:
                    non_synonymous += 1
                else:
                    synonymous += 1
        
        # Calculate dN/dS ratio
        dn_ds = non_synonymous / synonymous if synonymous > 0 else float('inf')
        
        # Save results
        results = {
            'non_synonymous': non_synonymous,
            'synonymous': synonymous,
            'dN/dS': dn_ds,
            'alignment_length': len(codon_alignment[0].seq)
        }
        
        return results
        
    except Exception as e:
        logging.error(f"Error calculating dN/dS: {str(e)}")
        raise

def run_dnds_analysis():
    """Run dN/dS analysis"""
    try:
        # Create output directory
        output_dir = Path('results/selection_analysis')
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Calculate dN/dS
        logging.info("Starting dN/dS calculation...")
        results = calculate_dn_ds('data/tp53_alignment.phy', 'results/tree/tree.nwk')
        
        # Save results
        with open(output_dir / 'dnds_results.json', 'w') as f:
            json.dump(results, f, indent=2)
        
        logging.info("dN/dS analysis completed successfully")
        return results
        
    except Exception as e:
        logging.error(f"Error in dN/dS analysis: {str(e)}")
        raise

def analyze_splicing_patterns(alignment_file):
    """Analyze splicing patterns using BioPython"""
    try:
        # Read alignment
        alignment = AlignIO.read(alignment_file, 'phylip-relaxed')
        
        # Get sequences
        human_seq = str(alignment[0].seq)
        pig_seq = str(alignment[1].seq)
        
        # Define splicing motifs
        splice_motifs = {
            'donor': ['GT'],
            'acceptor': ['AG']
        }
        
        # Find splicing sites
        human_sites = {
            'donor': [],
            'acceptor': []
        }
        pig_sites = {
            'donor': [],
            'acceptor': []
        }
        
        # Search for donor sites (GT)
        for i in range(len(human_seq) - 1):
            if human_seq[i:i+2] in splice_motifs['donor']:
                human_sites['donor'].append(i)
            if pig_seq[i:i+2] in splice_motifs['donor']:
                pig_sites['donor'].append(i)
        
        # Search for acceptor sites (AG)
        for i in range(len(human_seq) - 1):
            if human_seq[i:i+2] in splice_motifs['acceptor']:
                human_sites['acceptor'].append(i)
            if pig_seq[i:i+2] in splice_motifs['acceptor']:
                pig_sites['acceptor'].append(i)
        
        # Calculate differences
        donor_changes = len(set(human_sites['donor']) ^ set(pig_sites['donor']))
        acceptor_changes = len(set(human_sites['acceptor']) ^ set(pig_sites['acceptor']))
        
        # Save results
        results = {
            'human': human_sites,
            'pig': pig_sites,
            'changes': {
                'donor': donor_changes,
                'acceptor': acceptor_changes,
                'total': donor_changes + acceptor_changes
            }
        }
        
        return results
        
    except Exception as e:
        logging.error(f"Error analyzing splicing patterns: {str(e)}")
        raise

def run_splicing_analysis():
    """Run splicing pattern analysis"""
    try:
        # Create output directory
        output_dir = Path('results/splicing_analysis')
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Analyze splicing patterns
        logging.info("Starting splicing pattern analysis...")
        results = analyze_splicing_patterns('data/tp53_alignment.phy')
        
        # Save results
        with open(output_dir / 'splicing_results.json', 'w') as f:
            json.dump(results, f, indent=2)
        
        logging.info("Splicing pattern analysis completed successfully")
        return str(output_dir)
        
    except Exception as e:
        logging.error(f"Error in splicing analysis: {str(e)}")
        raise

def main():
    try:
        setup_logging()
        
        # Run dN/dS analysis
        logging.info("Starting dN/dS analysis...")
        dnds_results = run_dnds_analysis()
        
        # Run splicing analysis
        logging.info("Starting splicing analysis...")
        splicing_dir = run_splicing_analysis()
        
        # Save summary
        summary = {
            'dnds_results': dnds_results,
            'splicing_results': str(splicing_dir),
            'analysis_date': str(datetime.datetime.now())
        }
        
        with open('results/selection_analysis/analysis_summary.json', 'w') as f:
            json.dump(summary, f, indent=2)
        
        logging.info("Analysis completed successfully")
        
    except Exception as e:
        logging.error(f"Error in main analysis: {str(e)}")
        raise

if __name__ == "__main__":
    main()

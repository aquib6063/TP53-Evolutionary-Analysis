#!/usr/bin/env python3
"""
TP53 Evolutionary Analysis Pipeline
----------------------------------
Performs comprehensive analysis of TP53 CTD evolution across species.
"""
import os
import sys
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple

# Project directories
BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data"
RAW_DATA = DATA_DIR / "raw"
PROCESSED_DATA = DATA_DIR / "processed"
RESULTS_DIR = BASE_DIR / "results"
SCRIPTS_DIR = BASE_DIR / "scripts"

# Ensure directories exist
for d in [PROCESSED_DATA, RESULTS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# Reference coordinates (1-based)
HUMAN_CTD_START = 1071
HUMAN_CTD_END = 1179

class TP53Analyzer:
    def __init__(self):
        self.species_data = {}
        self.results = {}
        
    def run_command(self, cmd: str) -> str:
        """Execute shell command and return output."""
        try:
            result = subprocess.run(
                cmd, shell=True, check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            return result.stdout
        except subprocess.CalledProcessError as e:
            print(f"Error running command: {cmd}")
            print(f"Error: {e.stderr}")
            raise

    def analyze_composition(self):
        """Run base composition analysis."""
        from Bio import SeqIO
        from Bio.SeqUtils import GC, gc_fraction
        
        print("\n=== Running Base Composition Analysis ===")
        
        results = []
        for fasta in RAW_DATA.glob("*.fasta"):
            species = fasta.stem
            print(f"Processing {species}...")
            
            record = next(SeqIO.parse(fasta, "fasta"))
            cds = str(record.seq).upper()
            ctd = cds[HUMAN_CTD_START-1:HUMAN_CTD_END]  # Convert to 0-based
            
            # Calculate GC content
            cds_gc = GC(cds)
            ctd_gc = GC(ctd)
            delta_gc = ctd_gc - cds_gc
            
            # Calculate base composition
            def base_comp(seq):
                total = len(seq)
                return {b: seq.count(b)/total*100 for b in 'ATCG'}
            
            cds_comp = base_comp(cds)
            ctd_comp = base_comp(ctd)
            
            results.append({
                'species': species,
                'cds_length': len(cds),
                'cds_gc': cds_gc,
                'ctd_gc': ctd_gc,
                'delta_gc': delta_gc,
                **{f'cds_{k}': v for k, v in cds_comp.items()},
                **{f'ctd_{k}': v for k, v in ctd_comp.items()}
            })
        
        # Save results
        import pandas as pd
        df = pd.DataFrame(results)
        output_file = RESULTS_DIR / 'composition_analysis.csv'
        df.to_csv(output_file, index=False)
        print(f"\nComposition analysis saved to: {output_file}")
        return df

    def run_phylogeny(self):
        """Build phylogenetic tree using MAFFT and IQ-TREE."""
        print("\n=== Running Phylogenetic Analysis ===")
        
        # 1. Align sequences with MAFFT
        aligned_file = PROCESSED_DATA / 'tp53_aligned.fasta'
        if not aligned_file.exists():
            cmd = f"mafft --auto --thread -1 {RAW_DATA}/*.fasta > {aligned_file}"
            self.run_command(cmd)
        
        # 2. Build tree with IQ-TREE
        tree_file = RESULTS_DIR / 'tp53_phylogeny.tree'
        if not tree_file.exists():
            cmd = f"iqtree -s {aligned_file} -m MFP -bb 1000 -nt AUTO -pre {RESULTS_DIR/'tp53'}"
            self.run_command(cmd)
        
        print(f"Phylogenetic tree saved to: {tree_file}")
        return tree_file

    def predict_ptms(self):
        """Predict post-translational modifications."""
        print("\n=== Running PTM Prediction ===")
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        
        ptm_results = []
        
        for fasta in RAW_DATA.glob("*.fasta"):
            species = fasta.stem
            print(f"Predicting PTMs for {species}...")
            
            record = next(SeqIO.parse(fasta, "fasta"))
            cds = str(record.seq).upper()
            protein = str(Seq(cds).translate())
            
            # Simple PTM prediction (replace with actual PTM prediction tools)
            ptms = {
                'phosphorylation': len([i for i, aa in enumerate(protein) 
                                     if aa in 'STY' and i >= 100]),  # Example: S/T/Y after position 100
                'acetylation': len([i for i, aa in enumerate(protein) 
                                 if aa == 'K' and i >= 100]),
                'ubiquitination': len([i for i, aa in enumerate(protein) 
                                    if aa == 'K' and i >= 100])
            }
            
            ptm_results.append({
                'species': species,
                'protein_length': len(protein),
                **ptms
            })
        
        # Save results
        import pandas as pd
        df = pd.DataFrame(ptm_results)
        output_file = RESULTS_DIR / 'ptm_analysis.csv'
        df.to_csv(output_file, index=False)
        print(f"\nPTM analysis saved to: {output_file}")
        return df

def main():
    """Run the complete analysis pipeline."""
    analyzer = TP53Analyzer()
    
    # Run analyses
    comp_results = analyzer.analyze_composition()
    tree_file = analyzer.run_phylogeny()
    ptm_results = analyzer.predict_ptms()
    
    print("\n=== Analysis Complete ===")
    print(f"Results saved to: {RESULTS_DIR}")

if __name__ == "__main__":
    main()

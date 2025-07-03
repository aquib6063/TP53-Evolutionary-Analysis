import os
from Bio import SeqIO
import pandas as pd

# Base directory
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DIR = os.path.join(BASE_DIR, 'data', 'raw_sequences')
GC_DIR = os.path.join(BASE_DIR, 'results', 'gc_content')

def calculate_gc_content():
    """Calculate GC content for each species"""
    gc_results = []
    
    for fasta_file in os.listdir(RAW_DIR):
        if not fasta_file.endswith('.fasta'):
            continue
            
        species = fasta_file.replace('_TP53.fasta', '')
        
        with open(os.path.join(RAW_DIR, fasta_file)) as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                seq = str(record.seq)
                gc = (seq.count('G') + seq.count('C')) / len(seq) * 100
                gc_results.append({
                    'species': species,
                    'gc_content': gc
                })
    
    df = pd.DataFrame(gc_results)
    df.to_csv(os.path.join(GC_DIR, 'gc_by_domain.csv'), index=False)

if __name__ == "__main__":
    calculate_gc_content()

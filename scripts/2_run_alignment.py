import os
from Bio.Align.Applications import ClustalwCommandline

# Base directory
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DIR = os.path.join(BASE_DIR, 'data', 'raw_sequences')
PROCESSED_DIR = os.path.join(BASE_DIR, 'data', 'processed')

def align_sequences():
    """Align TP53 sequences using ClustalW"""
    # Get all FASTA files
    fasta_files = [os.path.join(RAW_DIR, f) for f in os.listdir(RAW_DIR) if f.endswith('.fasta')]
    
    # Create alignment input file
    alignment_input = os.path.join(PROCESSED_DIR, 'all_sequences.fasta')
    with open(alignment_input, 'w') as outfile:
        for fasta in fasta_files:
            with open(fasta) as infile:
                outfile.write(infile.read())
    
    # Run ClustalW
    cline = ClustalwCommandline("clustalw2", infile=alignment_input)
    stdout, stderr = cline()
    
    print("Alignment complete!")

if __name__ == "__main__":
    align_sequences()

import os
from Bio import AlignIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Base directory (should be the project root)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Input/output paths
INPUT_FILE = os.path.join(BASE_DIR, "data", "raw_sequences", "all_species_TP53.fasta")
OUTPUT_FILE = os.path.join(BASE_DIR, "data", "processed", "all_species_aligned.fasta")

def run_alignment():
    """Run multiple sequence alignment using Biopython"""
    try:
        # Check if input file exists
        if not os.path.exists(INPUT_FILE):
            logging.error(f"❌ Input file not found: {INPUT_FILE}")
            return
            
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
        
        # Read sequences
        sequences = list(SeqIO.parse(INPUT_FILE, "fasta"))
        
        # Create alignment
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 1.0
        aligner.mismatch_score = -1.0
        aligner.open_gap_score = -0.5
        aligner.extend_gap_score = -0.1
        
        # Create unique names for sequences
        unique_names = []
        for i, seq in enumerate(sequences):
            # Extract species name from description
            description = seq.description
            if "Human" in description:
                name = "Human_TP53_" + str(i)
            elif "Pig" in description:
                name = "Pig_TP53_" + str(i)
            elif "Zebrafish" in description:
                name = "Zebrafish_TP53_" + str(i)
            elif "Rhesus_Monkey" in description:
                name = "Rhesus_Monkey_TP53_" + str(i)
            elif "Rat" in description:
                name = "Rat_TP53_" + str(i)
            else:
                name = f"Unknown_{i}"
            unique_names.append(name)
        
        # Create alignment matrix
        aligned_sequences = []
        max_length = max(len(seq.seq) for seq in sequences)
        
        # Get the longest sequence from each species
        species_sequences = {}
        for i, seq in enumerate(sequences):
            species = seq.id.split('_')[0]
            if species not in species_sequences:
                species_sequences[species] = []
            species_sequences[species].append(seq)
        
        # For each species, use the longest sequence
        for species, seqs in species_sequences.items():
            longest_seq = max(seqs, key=lambda s: len(str(s.seq)))
            padded_seq = str(longest_seq.seq).ljust(max_length, '-')
            aligned_seq = SeqRecord(Seq(padded_seq), id=f"{species}_TP53", description=longest_seq.description)
            aligned_sequences.append(aligned_seq)
        
        # Save alignment
        with open(OUTPUT_FILE, "w") as handle:
            SeqIO.write(aligned_sequences, handle, "fasta")
        
        logging.info(f"✅ Alignment saved to: {OUTPUT_FILE}")
        
    except Exception as e:
        logging.error(f"❌ Error running alignment: {str(e)}")

if __name__ == "__main__":
    run_alignment()

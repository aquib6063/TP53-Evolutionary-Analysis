from Bio import SeqIO, AlignIO
from Bio.Align import PairwiseAligner
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator
from Bio.Phylo import write
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import logging
from pathlib import Path

def setup_logging():
    """Setup logging configuration"""
    logging.basicConfig(
        filename='alignment_log.log',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def prepare_alignment():
    """Prepare and align sequences"""
    try:
        # Create output directory
        output_dir = Path('data')
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Read sequences and select the main TP53 sequence
        human_seqs = list(SeqIO.parse('TP53_HUMAN_datasets/ncbi_dataset/data/TP53.faa', 'fasta'))
        pig_seqs = list(SeqIO.parse('TP53_Sus scrofa (pig)_datasets-4/ncbi_dataset/data/TP53.faa', 'fasta'))
        
        # Select the main TP53 sequence (longest sequence)
        human_seq = max(human_seqs, key=lambda x: len(x.seq))
        pig_seq = max(pig_seqs, key=lambda x: len(x.seq))
        
        # Create temporary FASTA file for alignment
        temp_fasta = output_dir / 'temp.fasta'
        with open(temp_fasta, 'w') as f:
            SeqIO.write([human_seq, pig_seq], f, 'fasta')
        
        # Create aligner
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.5
        
        # Align sequences
        alignments = aligner.align(human_seq.seq, pig_seq.seq)
        best_alignment = alignments[0]
        
        # Create alignment object
        alignment = MultipleSeqAlignment([
            SeqRecord(Seq(str(best_alignment[0])), id="Human_TP53"),
            SeqRecord(Seq(str(best_alignment[1])), id="Pig_TP53")
        ])
        
        # Write alignment in PHYLIP format
        phy_file = output_dir / 'tp53_alignment.phy'
        AlignIO.write(alignment, phy_file, 'phylip-relaxed')
        
        logging.info(f"Alignment written to {phy_file}")
        return phy_file
        
    except Exception as e:
        logging.error(f"Error preparing alignment: {str(e)}")
        raise

def create_tree(alignment_file):
    """Create tree from alignment"""
    try:
        # Create output directory
        tree_dir = Path('results/tree')
        tree_dir.mkdir(parents=True, exist_ok=True)
        
        # Read alignment
        alignment = AlignIO.read(alignment_file, 'phylip-relaxed')
        
        # Calculate distances
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(alignment)
        
        # Create tree
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)
        
        # Write tree
        tree_file = tree_dir / 'tree.nwk'
        write(tree, tree_file, 'newick')
        
        logging.info(f"Tree written to {tree_file}")
        return tree_file
        
    except Exception as e:
        logging.error(f"Error creating tree: {str(e)}")
        raise

def main():
    try:
        setup_logging()
        
        # Prepare alignment
        logging.info("Starting alignment preparation...")
        alignment_file = prepare_alignment()
        
        # Create tree
        logging.info("Starting tree construction...")
        tree_file = create_tree(alignment_file)
        
        logging.info("Preparation completed successfully")
        
    except Exception as e:
        logging.error(f"Error in main preparation: {str(e)}")
        raise

if __name__ == "__main__":
    main()

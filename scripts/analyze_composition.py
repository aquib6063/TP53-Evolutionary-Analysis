t#!/usr/bin/env python3
import os
import sys
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import pandas as pd
from pathlib import Path

# Mapping of sequence IDs to species (scientific name and common name)
SPECIES_MAP = {
    'NM_000546.5': {'scientific': 'Homo sapiens', 'common': 'Human'},
    'NM_001009294.1': {'scientific': 'Felis catus', 'common': 'Cat'},
    'NM_001047151.2': {'scientific': 'Macaca mulatta', 'common': 'Rhesus macaque'},
    'NM_001202405.1': {'scientific': 'Equus caballus', 'common': 'Horse'},
    'NM_001328588.1': {'scientific': 'Danio rerio', 'common': 'Zebrafish'},
    'NM_001389218.1': {'scientific': 'Canis lupus familiaris', 'common': 'Dog'},
    'NM_011640.3': {'scientific': 'Mus musculus', 'common': 'Mouse'},
    'NM_030989.4': {'scientific': 'Rattus norvegicus', 'common': 'Rat'},
    'NM_174201.2': {'scientific': 'Bos taurus', 'common': 'Cow'},
    'NM_213824.3': {'scientific': 'Sus scrofa', 'common': 'Pig'},
    'XM_003416902.4': {'scientific': 'Loxodonta africana', 'common': 'African elephant (predicted)'}
}

def get_species_info(seq_id):
    """Get species information from sequence ID."""
    return SPECIES_MAP.get(seq_id, {'scientific': 'Unknown', 'common': 'Unknown'})

def analyze_sequence(seq_record):
    """Calculate base composition statistics for a sequence."""
    seq = str(seq_record.seq).upper()
    total = len(seq)
    if total == 0:
        return {}
    
    # Get species information
    species_info = get_species_info(seq_record.id)
    
    # Calculate GC content using gc_fraction
    gc_content = gc_fraction(seq) * 100  # Convert to percentage
    
    # Calculate base counts
    counts = {
        'A': seq.count('A'),
        'T': seq.count('T'),
        'C': seq.count('C'),
        'G': seq.count('G'),
        'N': seq.count('N'),
        '-': seq.count('-')
    }
    
    # Calculate percentages
    total_bases = sum(counts.values())
    if total_bases == 0:
        return {}
    
    stats = {
        'Accession': seq_record.id,
        'Species_Scientific': species_info['scientific'],
        'Species_Common': species_info['common'],
        'Length': total_bases,
        'GC_percent': round(gc_content, 2),
        'A_percent': round((counts['A'] / total_bases) * 100, 2),
        'T_percent': round((counts['T'] / total_bases) * 100, 2),
        'C_percent': round((counts['C'] / total_bases) * 100, 2),
        'G_percent': round((counts['G'] / total_bases) * 100, 2),
        'N_or_gap_percent': round(((counts['N'] + counts['-']) / total_bases) * 100, 2)
    }
    
    return stats

def main():
    # Check if correct number of arguments are provided
    if len(sys.argv) != 3:
        print(f"Usage: python {sys.argv[0]} <input.fasta> <output.csv>")
        print(f"Current working directory: {os.getcwd()}")
        print("\nAvailable files in current directory:")
        os.system("ls -la")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Check if input file exists
    if not os.path.isfile(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        print(f"Current working directory: {os.getcwd()}")
        print("\nAvailable files in current directory:")
        os.system("ls -la")
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
            print(f"Created directory: {output_dir}")
        except OSError as e:
            print(f"Error creating directory {output_dir}: {e}")
            sys.exit(1)
    
    # Process sequences
    results = []
    try:
        for record in SeqIO.parse(input_file, "fasta"):
            stats = analyze_sequence(record)
            if stats:  # Only add if sequence was valid
                results.append(stats)
        
        # Convert to DataFrame and save
        if results:
            df = pd.DataFrame(results)
            # Reorder columns
            columns = ['Accession', 'Species_Scientific', 'Species_Common'] + \
                     [col for col in df.columns if col not in 
                      ['Accession', 'Species_Scientific', 'Species_Common']]
            df = df[columns]
            df.to_csv(output_file, index=False, float_format='%.2f')
            print(f"Analysis complete. Results saved to {output_file}")
            
            # Print summary statistics
            print("\nSummary Statistics by Species:")
            summary = df[['Species_Common', 'Length', 'GC_percent']].groupby('Species_Common').mean().round(2)
            print(summary)
            
        else:
            print("No valid sequences found in the input file.")
            
    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

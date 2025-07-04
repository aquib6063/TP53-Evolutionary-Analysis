import pandas as pd
import sys

def create_vcf(input_file, output_file):
    """Convert pig-specific variants to VCF format"""
    # Read variations
    variations = pd.read_csv(input_file)
    
    # Create VCF header
    header = """
##fileformat=VCFv4.2
##source=pig_specific_variants
##INFO=<ID=IMPACT,Number=1,Type=String,Description="Impact of the variant">
##INFO=<ID=TYPE,Number=1,Type=String,Description="Type of variant">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
"""
    
    # Convert variations to VCF format
    vcf_lines = []
    for _, row in variations.iterrows():
        # Extract position and bases
        pos = row['Position']
        ref = row['Human_Base']
        alt = row['Pig_Base']
        
        # Create VCF line
        vcf_line = f"chr17\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\tIMPACT=MODIFIER;TYPE=SNP"
        vcf_lines.append(vcf_line)
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write(header)
        f.write('\n'.join(vcf_lines))
    
    print(f"VCF file created: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python create_vcf.py <input_variants.txt> <output_vcf.vcf>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    create_vcf(input_file, output_file)

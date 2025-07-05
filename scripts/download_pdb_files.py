import requests
import logging
from pathlib import Path

# Setup logging
def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def download_pdb(pdb_id, output_file):
    """Download PDB file from RCSB"""
    try:
        url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
        response = requests.get(url)
        
        if response.status_code == 200:
            with open(output_file, 'w') as f:
                f.write(response.text)
            logging.info(f"Downloaded {pdb_id} to {output_file}")
        else:
            raise Exception(f"Failed to download {pdb_id}: HTTP {response.status_code}")
            
    except Exception as e:
        logging.error(f"Error downloading {pdb_id}: {str(e)}")
        raise

def main():
    try:
        setup_logging()
        
        # Create output directory
        output_dir = Path('data/structures')
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Download human TP53 structure (PDB ID: 1TSR)
        human_pdb = output_dir / 'human_tp53.pdb'
        download_pdb('1TSR', human_pdb)
        
        # Download pig TP53 structure (PDB ID: 6T21)
        pig_pdb = output_dir / 'pig_tp53.pdb'
        download_pdb('6T21', pig_pdb)
        
        logging.info("PDB files downloaded successfully")
        
    except Exception as e:
        logging.error(f"Error in main: {str(e)}")
        raise

if __name__ == "__main__":
    main()

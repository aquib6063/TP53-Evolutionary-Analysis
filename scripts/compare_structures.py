import requests
import logging
import os
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser
from io import StringIO
import json
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def setup_logging():
    """Setup logging configuration"""
    logging.basicConfig(
        filename='structure_analysis.log',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def fetch_alphafold_structure(uniprot_id: str) -> str:
    """Fetch AlphaFold structure from EBI"""
    try:
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
        response = requests.get(url)
        response.raise_for_status()
        return response.text
    except requests.RequestException as e:
        logging.error(f"Error fetching structure for {uniprot_id}: {str(e)}")
        raise

def compare_structures(human_pdb: str, pig_pdb: str, output_dir: Path):
    """Compare and visualize structures using matplotlib"""
    try:
        # Parse PDB structures
        parser = PDBParser(QUIET=True)
        
        # Create figure
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        # Process human structure
        human_structure = parser.get_structure('human', StringIO(human_pdb))
        human_atoms = []
        for atom in human_structure.get_atoms():
            if atom.get_name() == 'CA':  # Use alpha carbons
                human_atoms.append(atom.get_coord())
        human_atoms = np.array(human_atoms)
        
        # Process pig structure
        pig_structure = parser.get_structure('pig', StringIO(pig_pdb))
        pig_atoms = []
        for atom in pig_structure.get_atoms():
            if atom.get_name() == 'CA':  # Use alpha carbons
                pig_atoms.append(atom.get_coord())
        pig_atoms = np.array(pig_atoms)
        
        # Plot structures
        ax.plot(human_atoms[:, 0], human_atoms[:, 1], human_atoms[:, 2], 
                color='blue', label='Human TP53')
        ax.plot(pig_atoms[:, 0], pig_atoms[:, 1], pig_atoms[:, 2], 
                color='red', label='Pig TP53')
        
        # Set plot properties
        ax.set_title('TP53 Structure Comparison (Human vs Pig)')
        ax.legend()
        ax.set_xlabel('X coordinate')
        ax.set_ylabel('Y coordinate')
        ax.set_zlabel('Z coordinate')
        
        # Save figure
        img_path = output_dir / "structure_comparison.png"
        plt.savefig(img_path, dpi=300)
        plt.close()
        
        # Log basic structure information
        logging.info("Structures fetched successfully")
        
        return img_path
    except Exception as e:
        logging.error(f"Error in structure comparison: {str(e)}")
        raise

def analyze_local_structure_changes(human_pdb: str, pig_pdb: str, output_dir: Path):
    """Analyze local structural changes using matplotlib"""
    try:
        # Parse PDB structures
        parser = PDBParser(QUIET=True)
        
        # Create figure
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        # Process human structure
        human_structure = parser.get_structure('human', StringIO(human_pdb))
        human_atoms = []
        for atom in human_structure.get_atoms():
            if atom.get_name() == 'CA':  # Use alpha carbons
                human_atoms.append(atom.get_coord())
        human_atoms = np.array(human_atoms)
        
        # Process pig structure
        pig_structure = parser.get_structure('pig', StringIO(pig_pdb))
        pig_atoms = []
        for atom in pig_structure.get_atoms():
            if atom.get_name() == 'CA':  # Use alpha carbons
                pig_atoms.append(atom.get_coord())
        pig_atoms = np.array(pig_atoms)
        
        # Plot structures
        ax.plot(human_atoms[:, 0], human_atoms[:, 1], human_atoms[:, 2], 
                color='blue', label='Human TP53')
        ax.plot(pig_atoms[:, 0], pig_atoms[:, 1], pig_atoms[:, 2], 
                color='red', label='Pig TP53')
        
        # Highlight C-terminal domain (residues 360+ in human numbering)
        human_ctd_start = 360
        human_ctd_end = 400
        pig_ctd_start = 360
        pig_ctd_end = 400
        
        # Plot C-terminal domain
        ax.plot(human_atoms[human_ctd_start:human_ctd_end, 0], 
                human_atoms[human_ctd_start:human_ctd_end, 1],
                human_atoms[human_ctd_start:human_ctd_end, 2],
                color='blue', linewidth=3, label='Human CTD')
        
        ax.plot(pig_atoms[pig_ctd_start:pig_ctd_end, 0],
                pig_atoms[pig_ctd_start:pig_ctd_end, 1],
                pig_atoms[pig_ctd_start:pig_ctd_end, 2],
                color='red', linewidth=3, label='Pig CTD')
        
        # Set plot properties
        ax.set_title('TP53 C-terminal Domain Comparison')
        ax.legend()
        ax.set_xlabel('X coordinate')
        ax.set_ylabel('Y coordinate')
        ax.set_zlabel('Z coordinate')
        
        # Save figure
        img_path = output_dir / "ctd_structure_comparison.png"
        plt.savefig(img_path, dpi=300)
        plt.close()
        
        return img_path
    except Exception as e:
        logging.error(f"Error in local structure analysis: {str(e)}")
        raise

def main():
    try:
        # Setup
        setup_logging()
        output_dir = Path(__file__).parent.parent / "results" / "structure_analysis"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Fetch structures
        logging.info("Fetching AlphaFold structures...")
        human_pdb = fetch_alphafold_structure("P04637")
        pig_pdb = fetch_alphafold_structure("A0A287BKT5")
        
        # Compare full structures
        logging.info("Comparing full structures...")
        full_comparison = compare_structures(human_pdb, pig_pdb, output_dir)
        logging.info(f"Saved full structure comparison to {full_comparison}")
        
        # Analyze C-terminal domain
        logging.info("Analyzing C-terminal domain...")
        ctd_comparison = analyze_local_structure_changes(human_pdb, pig_pdb, output_dir)
        logging.info(f"Saved C-terminal domain comparison to {ctd_comparison}")
        
        # Save results
        results = {
            'full_comparison': str(full_comparison),
            'ctd_comparison': str(ctd_comparison)
        }
        
        with open(output_dir / 'structure_analysis_results.json', 'w') as f:
            json.dump(results, f, indent=2)
        
        logging.info("Structure analysis completed successfully")
        
    except Exception as e:
        logging.error(f"Error in main analysis: {str(e)}")
        raise

if __name__ == "__main__":
    main()

import os
import logging
from pathlib import Path

# Setup logging
def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def write_pymol_script():
    """Write PyMOL script for structure comparison"""
    try:
        script_content = """
# Load structures
load human_tp53.pdb
load pig_tp53.pdb

# Align structures
align pig_tp53, human_tp53

# Visualize
spectrum b, rainbow, human_tp53
show cartoon, human_tp53
show surface, pig_tp53, transparency=0.6

# Highlight CTD
select human_ctd, resi 360-393
select pig_ctd, resi 360-393
show sticks, human_ctd
show mesh, pig_ctd
color red, human_ctd
color blue, pig_ctd

# Save image
ray 1600, 1200
png structural_comparison.png
"""

        # Create output directory
        output_dir = Path('results/figures')
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Write script
        script_path = output_dir / 'structure_comparison.pml'
        with open(script_path, 'w') as f:
            f.write(script_content)
            
        logging.info(f"PyMOL script written to: {script_path}")
        return script_path
        
    except Exception as e:
        logging.error(f"Error writing PyMOL script: {str(e)}")
        raise

def run_pymol(script_path):
    """Run PyMOL with the script"""
    try:
        # Run PyMOL
        os.system(f'pymol -c {script_path}')
        
        # Move the generated image to results directory
        output_dir = Path('results/figures')
        image_path = output_dir / 'structural_comparison.png'
        if image_path.exists():
            logging.info(f"Image generated: {image_path}")
        else:
            logging.warning("No image was generated")
            
    except Exception as e:
        logging.error(f"Error running PyMOL: {str(e)}")
        raise

def main():
    try:
        setup_logging()
        logging.info("Starting structure visualization...")
        
        # Write and run PyMOL script
        script_path = write_pymol_script()
        run_pymol(script_path)
        
        logging.info("Visualization completed successfully")
        
    except Exception as e:
        logging.error(f"Error in main: {str(e)}")
        raise

if __name__ == "__main__":
    main()

import logging
import os
from Bio import SeqIO
from Bio.Phylo.PAML import codeml
from pathlib import Path
import pandas as pd
import numpy as np

# Setup logging
def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def create_codeml_config():
    """Create CodeML configuration files"""
    try:
        # Create output directory
        output_dir = Path('results/selection_analysis')
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create M0 (one-ratio) configuration
        m0_config = {
            'seqfile': 'data/sequences/aligned_codons.fasta',
            'treefile': 'results/alignment_and_tree/tree.nwk',
            'outfile': 'results/selection_analysis/m0_results.txt',
            'noisy': 0,
            'verbose': 1,
            'runmode': 0,
            'seqtype': 1,
            'CodonFreq': 2,
            'clock': 0,
            'aaDist': 0,
            'model': 0,
            'NSsites': [0],
            'icode': 0,
            'fix_kappa': 0,
            'kappa': 2,
            'fix_omega': 0,
            'omega': 0.4,
            'fix_alpha': 1,
            'alpha': 0,
            'Malpha': 0,
            'ncatG': 10,
            'getSE': 0,
            'RateAncestor': 0,
            'Small_Diff': 0.5e-6,
            'cleandata': 1
        }
        
        # Create M8 (positive selection) configuration
        m8_config = m0_config.copy()
        m8_config.update({
            'NSsites': [8],
            'fix_omega': 0,
            'omega': 1.5,
            'fix_alpha': 1,
            'alpha': 0,
            'Malpha': 0,
            'ncatG': 10
        })
        
        return m0_config, m8_config
        
    except Exception as e:
        logging.error(f"Error creating CodeML config: {str(e)}")
        raise

def run_codeml(config, model_name):
    """Run CodeML with given configuration"""
    try:
        # Create output directory
        output_dir = Path('results/selection_analysis')
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create control file
        ctl_file = output_dir / f'{model_name}.ctl'
        with open(ctl_file, 'w') as f:
            for key, value in config.items():
                if isinstance(value, list):
                    f.write(f'{key} = {" ".join(map(str, value))}\n')
                else:
                    f.write(f'{key} = {value}\n')
        
        # Run CodeML
        cml = codeml.Codeml()
        cml.set_options(**config)
        results = cml.run()
        
        # Save results
        results_file = output_dir / f'{model_name}_results.txt'
        with open(results_file, 'w') as f:
            f.write(str(results))
        
        logging.info(f"CodeML {model_name} results saved to: {results_file}")
        return results
        
    except Exception as e:
        logging.error(f"Error running CodeML {model_name}: {str(e)}")
        raise

def perform_likelihood_ratio_test(m0_results, m8_results):
    """Perform Likelihood Ratio Test between M0 and M8 models"""
    try:
        # Extract log-likelihood values
        lnL_m0 = m0_results.get('NSsites')
        lnL_m8 = m8_results.get('NSsites')
        
        if lnL_m0 is None or lnL_m8 is None:
            raise ValueError("Log-likelihood values not found in results")
            
        # Calculate LRT statistic
        LRT = -2 * (lnL_m0 - lnL_m8)
        
        # Calculate p-value (2 degrees of freedom)
        p_value = 1 - stats.chi2.cdf(LRT, df=2)
        
        # Determine significance
        significant = p_value < 0.05
        
        return {
            'LRT_statistic': LRT,
            'p_value': p_value,
            'significant': significant,
            'positive_selection': m8_results.get('site_classes', {}).get('2', 0)
        }
        
    except Exception as e:
        logging.error(f"Error performing LRT: {str(e)}")
        raise

def analyze_domain_selection():
    """Analyze selection pressure in each domain"""
    try:
        # Create output directory
        output_dir = Path('results/selection_analysis')
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Domain boundaries
        domain_boundaries = {
            'N-terminal': (0, 42*3),
            'DNA-binding': (94*3, 289*3),
            'Oligomerization': (323*3, 355*3),
            'C-terminal': (363*3, 393*3)
        }
        
        # Process each domain
        results = []
        for domain, (start, end) in domain_boundaries.items():
            # Extract domain sequences
            sequences = []
            for record in SeqIO.parse('data/sequences/aligned_codons.fasta', 'fasta'):
                domain_seq = str(record.seq)[start:end]
                sequences.append(f'>{record.id}\n{domain_seq}')
            
            # Write domain sequences
            domain_file = output_dir / f'{domain}_sequences.fasta'
            with open(domain_file, 'w') as f:
                f.write('\n'.join(sequences))
            
            # Update config for domain analysis
            m0_config, m8_config = create_codeml_config()
            m0_config['seqfile'] = str(domain_file)
            m8_config['seqfile'] = str(domain_file)
            
            # Run analysis
            m0_results = run_codeml(m0_config, f'm0_{domain}')
            m8_results = run_codeml(m8_config, f'm8_{domain}')
            
            # Perform LRT
            lrt_results = perform_likelihood_ratio_test(m0_results, m8_results)
            
            # Save results
            results.append({
                'domain': domain,
                **lrt_results
            })
        
        # Save domain-specific results
        results_df = pd.DataFrame(results)
        output_file = output_dir / 'domain_selection_analysis.csv'
        results_df.to_csv(output_file, index=False)
        logging.info(f"Domain selection analysis saved to: {output_file}")
        
        return results_df
        
    except Exception as e:
        logging.error(f"Error analyzing domain selection: {str(e)}")
        raise

def main():
    try:
        setup_logging()
        logging.info("Starting selection pressure analysis...")
        
        # Run full sequence analysis
        m0_config, m8_config = create_codeml_config()
        m0_results = run_codeml(m0_config, 'm0_full')
        m8_results = run_codeml(m8_config, 'm8_full')
        
        # Perform LRT
        lrt_results = perform_likelihood_ratio_test(m0_results, m8_results)
        logging.info(f"LRT Results: {lrt_results}")
        
        # Analyze domains
        domain_results = analyze_domain_selection()
        
        logging.info("Selection pressure analysis completed successfully")
        
    except Exception as e:
        logging.error(f"Error in main: {str(e)}")
        raise

if __name__ == "__main__":
    main()

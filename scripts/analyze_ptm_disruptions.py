import pandas as pd
from Bio import ExPASy, SeqIO, Align
import requests
import logging
import os
from typing import Dict, List, Optional
from collections import defaultdict
# NetPhos integration is currently not available
# from bioservices import NetPhos

# Define PTM residues and their possible modifications
PTM_RESIDUES = {
    'S': ['phosphorylation'],
    'T': ['phosphorylation'],
    'Y': ['phosphorylation'],
    'K': ['acetylation', 'methylation', 'ubiquitination'],
    'R': ['methylation'],
    'E': ['ubiquitination'],
    'N': ['ubiquitination'],
    'C': ['ubiquitination']
}

# Define PTM prediction thresholds (more sensitive)
PTM_THRESHOLDS = {
    'phosphorylation': 0.5,  # Lowered from 0.8
    'acetylation': 0.4,     # Lowered from 0.7
    'methylation': 0.3,     # Lowered from 0.6
    'ubiquitination': 0.2   # Lowered from 0.5
}

# Configure logging
logging.basicConfig(
    filename='ptm_analysis.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Base directory
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Output directory
OUTPUT_DIR = os.path.join(BASE_DIR, "results", "ptm_analysis")

def get_sequence(uniprot_id: str):
    """Fetch sequence from UniProt"""
    try:
        handle = ExPASy.get_sprot_raw(uniprot_id)
        sequence = str(SeqIO.read(handle, "swiss").seq)
        logging.info(f"Successfully fetched sequence for {uniprot_id}")
        return sequence
    except Exception as e:
        logging.error(f"Error fetching sequence for {uniprot_id}: {str(e)}")
        raise

def fetch_ptm_sites(uniprot_id: str) -> Dict[int, Dict]:
    """Fetch PTM sites from multiple sources"""
    try:
        ptm_sites = {}
        
        # Try UniProt API
        try:
            # UniProt API endpoint
            url = f"https://www.uniprot.org/uniprotkb/{uniprot_id}.xml"
            
            # Get the XML content
            response = requests.get(url)
            response.raise_for_status()
            
            # Parse the XML to extract PTM sites
            from xml.etree import ElementTree as ET
            root = ET.fromstring(response.text)
            
            # Find all features (PTMs)
            for feature in root.findall('.//feature'):
                try:
                    # Get feature type
                    feature_type = feature.get('type')
                    
                    # Get position
                    position = int(feature.find('position').get('position'))
                    
                    # Get residue
                    residue = feature.find('original').text
                    
                    # Get modification type
                    description = feature.find('description').text
                    
                    # Map description to PTM type
                    ptm_type = 'unknown'
                    if 'phosphorylation' in description.lower():
                        ptm_type = 'phosphorylation'
                    elif 'acetylation' in description.lower():
                        ptm_type = 'acetylation'
                    elif 'ubiquitination' in description.lower():
                        ptm_type = 'ubiquitination'
                    elif 'methylation' in description.lower():
                        ptm_type = 'methylation'
                    elif 'glycosylation' in description.lower():
                        ptm_type = 'glycosylation'
                    
                    ptm_sites[position] = {
                        'residue': residue,
                        'site_type': ptm_type,
                        'description': description,
                        'source': 'uniprot'
                    }
                    logging.info(f"[UniProt] Found PTM site at position {position}: {ptm_type} - {description}")
                except Exception as e:
                    logging.warning(f"Error parsing UniProt PTM site: {str(e)}")
                    continue
        except Exception as e:
            logging.warning(f"Error fetching from UniProt: {str(e)}")
        
        # Try PhosphoSitePlus API
        try:
            # PhosphoSitePlus API endpoint
            url = f"https://www.phosphosite.org/uniprotSearchAction.action?uniprotAccession={uniprot_id}&siteType=all"
            
            # Get the page content
            response = requests.get(url)
            response.raise_for_status()
            
            # Parse the HTML to extract PTM sites
            from bs4 import BeautifulSoup
            soup = BeautifulSoup(response.text, 'html.parser')
            
            # Find all PTM sites
            site_elements = soup.find_all('div', class_='site')
            for element in site_elements:
                try:
                    # Extract position and residue
                    position = int(element.find('span', class_='position').text)
                    residue = element.find('span', class_='residue').text
                    
                    # Extract site type
                    site_type = element.find('span', class_='siteType').text
                    
                    # Add to PTM sites
                    ptm_sites[position] = {
                        'residue': residue,
                        'site_type': site_type,
                        'description': f"PhosphoSitePlus: {site_type}",
                        'source': 'phosphosite'
                    }
                    logging.info(f"[PhosphoSite] Found PTM site at position {position}: {site_type}")
                except Exception as e:
                    logging.warning(f"Error parsing PhosphoSite PTM site: {str(e)}")
                    continue
        except Exception as e:
            logging.warning(f"Error fetching from PhosphoSitePlus: {str(e)}")
        
        logging.info(f"Total PTM sites found for {uniprot_id}: {len(ptm_sites)}")
        return ptm_sites
    except Exception as e:
        logging.error(f"Error fetching PTM sites for {uniprot_id}: {str(e)}")
        raise

def local_ptm_prediction(sequence: str) -> Dict[int, float]:
    """Local PTM prediction using sequence context"""
    try:
        predictions = {}
        for i, res in enumerate(sequence):
            if res in PTM_RESIDUES:
                # Get local context (5 residues on each side)
                context_start = max(0, i-5)
                context_end = min(len(sequence), i+6)
                context = sequence[context_start:context_end]
                
                # Calculate context score
                context_score = sum(1 for r in context if r == res) / len(context)
                
                # Use more sensitive threshold
                if context_score > 0.2:  # Lowered from 0.3
                    predictions[i] = context_score
        
        return predictions
    except Exception as e:
        logging.error(f"Error in local PTM prediction: {str(e)}")
        return {}

def predict_ptm_sites(sequence, start, end, accession):
    """Predict PTM sites using local sequence context"""
    try:
        # Get local PTM predictions
        local_preds = local_ptm_prediction(sequence)
        
        # Initialize PTM predictions
        predictions = []
        
        # Analyze each residue
        for i in range(len(sequence)):
            res = sequence[i]
            pos = start + i
            
            # Skip if not a PTM residue
            if res not in PTM_RESIDUES:
                continue
                
            # Get local context (5 residues on each side)
            context_start = max(0, i-5)
            context_end = min(len(sequence), i+6)
            context = sequence[context_start:context_end]
            
            # Calculate conservation score based on context
            conservation = sum(1 for r in context if r == res) / len(context)
            
            # Get possible PTM types
            ptm_types = PTM_RESIDUES[res]
            
            # Get local prediction score if available
            local_score = local_preds.get(i, 0)
            
            # Calculate scores for each PTM type
            for ptm_type in ptm_types:
                threshold = PTM_THRESHOLDS[ptm_type]
                
                # Calculate context score based on amino acid properties
                context_score = 0
                for r in context:
                    if r == res:
                        context_score += 0.5  # Same residue
                    elif r in PTM_RESIDUES:
                        context_score += 0.3  # Other PTM residue
                    else:
                        context_score += 0.1  # Non-PTM residue
                
                context_score /= len(context)
                
                # Combine scores with local prediction
                if ptm_type == 'phosphorylation':
                    score = (conservation + context_score + local_score) / 3
                else:
                    score = (conservation + context_score) / 2
                
                # Add prediction if above threshold
                if score >= threshold:
                    predictions.append({
                        'position': pos,
                        'residue': res,
                        'ptm_type': ptm_type,
                        'score': score,
                        'conservation': conservation,
                        'context': context,
                        'method': 'local_context',
                        'local_score': local_score if ptm_type == 'phosphorylation' else None
                    })
        
        return predictions
    except Exception as e:
        logging.error(f"Error predicting PTM sites: {str(e)}")
        return []

def get_possible_ptms(residue):
    """Get possible PTM types for a given residue"""
    return PTM_RESIDUES.get(residue, [])

def analyze_ptm_disruptions(human_seq, pig_seq, start, end):
    """Analyze PTM disruptions between human and pig sequences"""
    try:
        logging.info(f"Starting PTM analysis for domain {start}-{end}")
        
        # Validate sequences
        if not human_seq or not pig_seq:
            raise ValueError(f"Sequence has zero length (human: {len(human_seq)}, pig: {len(pig_seq)})")
        
        # Get sequence similarity
        aligner = Align.PairwiseAligner()
        alignment = aligner.align(human_seq, pig_seq)
        score = alignment[0].score
        logging.info(f"Sequence similarity score: {score}")
        
        # Create DataFrame to store results
        results = []
        
        # Predict PTM sites for both sequences
        human_ptm_sites = predict_ptm_sites(human_seq, start, end, "P04637")
        pig_ptm_sites = predict_ptm_sites(pig_seq, start, end, "A0A287BKT5")
        
        logging.info(f"Human PTM sites found: {len(human_ptm_sites)}")
        logging.info(f"Pig PTM sites found: {len(pig_ptm_sites)}")
        
        # Convert predictions to dictionaries by position
        human_ptm_dict = {site['position']: site for site in human_ptm_sites}
        pig_ptm_dict = {site['position']: site for site in pig_ptm_sites}
        
        # Compare PTM sites
        for i in range(len(human_seq)):
            human_res = human_seq[i]
            pig_res = pig_seq[i]
            pos = start + i
            
            human_site = human_ptm_dict.get(pos, None)
            pig_site = pig_ptm_dict.get(pos, None)
            
            if human_site or pig_site:  # If either sequence has a predicted PTM
                results.append({
                    'position': pos,
                    'human_residue': human_res,
                    'pig_residue': pig_res,
                    'human_ptm_type': human_site['ptm_type'] if human_site else None,
                    'pig_ptm_type': pig_site['ptm_type'] if pig_site else None,
                    'human_score': human_site['score'] if human_site else 0,
                    'pig_score': pig_site['score'] if pig_site else 0,
                    'human_local_score': human_site['local_score'] if human_site and 'local_score' in human_site else None,
                    'pig_local_score': pig_site['local_score'] if pig_site and 'local_score' in pig_site else None,
                    'sequence_similarity': score,
                    'disruption_type': 'loss' if human_site and not pig_site else 
                                     'gain' if not human_site and pig_site else 
                                     'change' if human_site and pig_site and human_site['ptm_type'] != pig_site['ptm_type'] else None
                })
        
        return pd.DataFrame(results)
    except Exception as e:
        logging.error(f"Error analyzing PTM disruptions: {str(e)}")
        return pd.DataFrame()
    except Exception as e:
        logging.error(f"Error analyzing PTM disruptions: {str(e)}")
        raise

def main():
    """Main function to run the analysis"""
    try:
        # Set up logging
        logging.basicConfig(
            filename='ptm_analysis.log',
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        
        # Fetch sequences
        logging.info("Fetching sequences...")
        human_seq = get_sequence("P04637")  # Human TP53
        pig_seq = get_sequence("A0A287BKT5")  # Pig TP53
        
        # Define multiple domain boundaries to test
        domain_boundaries = [
            (360, 393),  # Standard CTD
            (350, 400),  # Extended CTD
            (325, 393),  # Including regulatory region
            (300, 450),  # Larger region
            (330, 380),  # Narrower region
            (340, 410),  # Alternative CTD
        ]
        
        all_results = []
        
        # Test each domain definition
        for i, (start, end) in enumerate(domain_boundaries):
            logging.info(f"\nTesting domain {i+1}: {start}-{end}")
            
            # Get domain sequences
            human_domain = human_seq[start:end]
            pig_domain = pig_seq[start:end]
            
            # Log domain sequences
            logging.info(f"Human domain sequence: {human_domain}")
            logging.info(f"Pig domain sequence: {pig_domain}")
            
            # Get sequence similarity
            aligner = Align.PairwiseAligner()
            alignment = aligner.align(human_domain, pig_domain)
            score = alignment[0].score
            logging.info(f"Sequence similarity score: {score}")
            
            # Predict PTM sites for both sequences
            logging.info(f"Predicting PTM sites for domain {start}-{end}...")
            human_ptm_sites = predict_ptm_sites(human_domain, start, end, "P04637")
            pig_ptm_sites = predict_ptm_sites(pig_domain, start, end, "A0A287BKT5")
            
            logging.info(f"Human PTM sites found: {len(human_ptm_sites)}")
            logging.info(f"Pig PTM sites found: {len(pig_ptm_sites)}")
            
            # Analyze PTM disruptions for this domain
            logging.info(f"Analyzing PTM disruptions for domain {start}-{end}...")
            disruptions = analyze_ptm_disruptions(human_domain, pig_domain, start, end)
            
            # Add domain information to results
            if not disruptions.empty:
                disruptions['domain_start'] = start
                disruptions['domain_end'] = end
                disruptions['domain_size'] = end - start
                disruptions['sequence_similarity'] = score
            
            all_results.append(disruptions)
            
            # Save results for this domain
            output_dir = "results/ptm_analysis"
            os.makedirs(output_dir, exist_ok=True)
            output_file = os.path.join(output_dir, f"ptm_disruptions_domain_{start}_{end}.csv")
            disruptions.to_csv(output_file, index=False)
            logging.info(f"Saved PTM disruption analysis for domain {start}-{end} to {output_file}")
            
        # Combine all results
        combined_results = pd.concat(all_results, ignore_index=True)
        combined_file = os.path.join(output_dir, "ptm_disruptions_combined.csv")
        combined_results.to_csv(combined_file, index=False)
        logging.info(f"Saved combined PTM disruption analysis to {combined_file}")
        
        # Generate summary statistics
        logging.info("\nPTM Analysis Summary:")
        logging.info(f"Total disrupted sites across all domains: {len(combined_results)}")
        
        if not combined_results.empty:
            # Analyze domain-specific patterns
            domain_stats = combined_results.groupby(['domain_start', 'domain_end']).agg({
                'human_ptm_type': lambda x: x.mode()[0] if len(x) > 0 else None,
                'pig_ptm_type': lambda x: x.mode()[0] if len(x) > 0 else None,
                'human_residue': lambda x: x.mode()[0] if len(x) > 0 else None,
                'pig_residue': lambda x: x.mode()[0] if len(x) > 0 else None,
                'sequence_similarity': 'mean'
            }).reset_index()
            
            print("\nDomain-Specific Statistics:")
            for _, row in domain_stats.iterrows():
                print(f"Domain {row['domain_start']}-{row['domain_end']}:")
                print("  Total disruptions:", len(combined_results[(combined_results['domain_start'] == row['domain_start']) & 
                                                   (combined_results['domain_end'] == row['domain_end'])]))
                print(f"  Average similarity: {row['sequence_similarity']:.2f}")
                print(f"  Human PTM type: {row['human_ptm_type']}")
                print(f"  Pig PTM type: {row['pig_ptm_type']}")
                print(f"  Human residue: {row['human_residue']}")
                print(f"  Pig residue: {row['pig_residue']}")
    except Exception as e:
        logging.error(f"Error in main analysis: {str(e)}")
        raise

if __name__ == "__main__":
    main()

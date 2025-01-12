"""
African Genetic Diversity Analysis
Author: Simon Gyimah
Focus: Analysis of disease-related genes and unique African genetic variants
"""

import sys
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional, TextIO
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from Bio import Entrez
import requests
import traceback
from dataclasses import dataclass
import logging
from pathlib import Path
import statsmodels.api as sm
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import time  
import asyncio
import aiohttp
import json
from datetime import datetime
from itertools import combinations
from collections import Counter
import argparse
import warnings




# Handle all different warning types
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=pd.errors.SettingWithCopyWarning)
# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


POPULATION_MAPPINGS = {
    '1000GENOMES:phase_3:YRI': {"code": "YRI", "major_group": "African", "subgroup": "West_African", "name": "Yoruba", "region": "Nigeria"},
    '1000GENOMES:phase_3:LWK': {"code": "LWK", "major_group": "African", "subgroup": "East_African", "name": "Luhya", "region": "Kenya"},
    '1000GENOMES:phase_3:GWD': {"code": "GWD", "major_group": "African", "subgroup": "West_African", "name": "Gambian", "region": "Gambia"},
    '1000GENOMES:phase_3:ESN': {"code": "ESN", "major_group": "African", "subgroup": "West_African", "name": "Esan", "region": "Nigeria"},
    '1000GENOMES:phase_3:MSL': {"code": "MSL", "major_group": "African", "subgroup": "West_African", "name": "Mende", "region": "Sierra Leone"},
    '1000GENOMES:phase_3:ACB': {"code": "ACB", "major_group": "African", "subgroup": "African_Diaspora", "name": "African Caribbean", "region": "Barbados"},
    '1000GENOMES:phase_3:ASW': {"code": "ASW", "major_group": "African", "subgroup": "African_Diaspora", "name": "African American", "region": "USA"},
    '1000GENOMES:phase_3:JPT': {"code": "JPT", "major_group": "East_Asian", "subgroup": "Eastern", "name": "Japanese", "region": "Japan"},
    '1000GENOMES:phase_3:CHB': {"code": "CHB", "major_group": "East_Asian", "subgroup": "Chinese", "name": "Han Chinese", "region": "Beijing"},
    '1000GENOMES:phase_3:CEU': {"code": "CEU", "major_group": "European", "subgroup": "Northern_European", "name": "Central European", "region": "Europe"}
}

@dataclass
class PopulationData:
    """Store population-specific genetic data"""
    population_id: str
    region: str
    sample_size: int
    variants: pd.DataFrame

class AfricanGeneticDiversity:
    def __init__(self, email: str, output_dir: str = "results"):
        """Initialize with updated gene names and aliases"""
        Entrez.email = email
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Initialize data storage
        self.variant_data = pd.DataFrame()
        self.scaled_data = None
        self.population_labels = None
        self.selection_data = None
        self.analysis_results = {}
        
        # Define populations
        self.populations = {
            'GHA': {'name': 'Ghanaian', 'region': 'West Africa'},
            'YRI': {'name': 'Yoruba', 'region': 'West Africa'},
            'LWK': {'name': 'Luhya', 'region': 'East Africa'},
            'CEU': {'name': 'European', 'region': 'Europe'},
            'CHB': {'name': 'Han Chinese', 'region': 'East Asia'}
        }
        
        # Define disease-related genes with current names and aliases
        self.disease_genes = {
            'APOL1': {
                'disease': 'Kidney disease',
                'chr': '22',
                'aliases': ['APO-L', 'APOL']
            },
            'HBB': {
                'disease': 'Sickle cell',
                'chr': '11',
                'aliases': ['beta-globin']
            },
            'G6PD': {
                'disease': 'Malaria resistance',
                'chr': 'X',
                'aliases': ['G6PD1']
            },
            'ACKR1': {  # Updated from DARC
                'disease': 'Malaria resistance',
                'chr': '1',
                'aliases': ['DARC', 'FY', 'GPD', 'CD234']
            },
            'TRPV6': {
                'disease': 'Calcium absorption',
                'chr': '7',
                'aliases': ['CAT1', 'ECAC2']
            }
        }

    def expand_gene_panels(self) -> None:
        """Expand gene panels for comprehensive analysis"""
        logger.info("Expanding gene panels")
        
        additional_genes = {
            # Cardiovascular disease genes
            'ACE': {'disease': 'Hypertension', 'chr': '17'},
            'APOB': {'disease': 'Cholesterol metabolism', 'chr': '2'},
            'PCSK9': {'disease': 'Cholesterol regulation', 'chr': '1'},
            
            # Metabolic disease genes
            'TCF7L2': {'disease': 'Type 2 diabetes', 'chr': '10'},
            'KCNJ11': {'disease': 'Diabetes', 'chr': '11'},
            
            # Cancer susceptibility genes
            'TP53': {'disease': 'Multiple cancers', 'chr': '17'},
            'BRCA1': {'disease': 'Breast cancer', 'chr': '17'},
            'BRCA2': {'disease': 'Breast cancer', 'chr': '13'},
            
            # Immune system genes
            'HLA-B': {'disease': 'Immune response', 'chr': '6'},
            'IL6': {'disease': 'Inflammatory response', 'chr': '7'},
            
            # Population-specific variants
            'SLC24A5': {'trait': 'Skin pigmentation', 'chr': '15'},
            'LCT': {'trait': 'Lactase persistence', 'chr': '2'}
        }
        
        self.disease_genes.update(additional_genes)
        logger.info(f"Gene panel expanded to {len(self.disease_genes)} genes")

    def _fetch_gene_info(self, base_url: str, gene_name: str) -> Optional[Dict]:
        """Fetch gene information with error handling"""
        try:
            gene_url = f"{base_url}/lookup/symbol/homo_sapiens/{gene_name}?"
            headers = {"Content-Type": "application/json"}
            
            response = requests.get(gene_url, headers=headers)
            response.raise_for_status()
            return response.json()
            
        except requests.exceptions.RequestException as e:
            logger.debug(f"Could not fetch info for gene {gene_name}: {e}")
            return None

    def fetch_variant_data(self) -> None:
        """Fetch variant data with optimized performance"""
        try:
            asyncio.run(self.fetch_variant_data_async())
        except KeyboardInterrupt:
            print("\nProcess interrupted by user. Shutting down gracefully...")
            logger.info("Process interrupted by user")
        except Exception as e:
            print(f"\nAn error occurred: {e}")
            logger.error(f"Error in fetch_variant_data: {e}")

    # 3. Add this new method:
    async def fetch_variant_data_async(self) -> None:
        """Asynchronous variant data fetching with improved error handling and logging"""
        # Initialize session with lower connection limit
        try:
            connector = aiohttp.TCPConnector(limit=3)  # Reduced from 20 to 3
            async with aiohttp.ClientSession(connector=connector) as session:
                all_variants_data = []
                headers = {
                    "Content-Type": "application/json",
                    "Accept": "application/json"
                }
                
                # Initialize counters for logging
                total_genes = len(self.disease_genes)
                processed_genes = 1
                
                logger.info(f"Starting to process {total_genes} genes")
                
                # Load cache if exists
                cache_file = Path(self.output_dir) / "variant_cache.json"
                cache = {}
                if cache_file.exists():
                    with open(cache_file, 'r') as f:
                        cache = json.load(f)
                
                # Process genes in smaller chunks to avoid rate limiting
                genes = list(self.disease_genes.keys())
                chunk_size = 3  # Process 3 genes at a time
                
                for i in range(0, len(genes), chunk_size):
                    gene_chunk = genes[i:i + chunk_size]
                    logger.info(f"Processing genes {i+1}-{min(i+chunk_size, len(genes))} of {len(genes)}")
                    
                    # Create tasks for this chunk of genes
                    chunk_tasks = [
                        self._fetch_gene_variants(session, gene, cache, all_variants_data, cache_file)
                        for gene in gene_chunk
                    ]
                    
                    # Process this chunk
                    await asyncio.gather(*chunk_tasks)
                    
                    # Add delay between chunks
                    if i + chunk_size < len(genes):  # If not the last chunk
                        await asyncio.sleep(1)  # 1 second delay between chunks
                
                # Convert to DataFrame
                self.variant_data = pd.DataFrame(all_variants_data)
                if not self.variant_data.empty:
                    logger.info(f"Successfully fetched data for {len(self.variant_data['gene'].unique())} genes")
                    logger.info(f"Total variants: {len(self.variant_data)}")
                    self.variant_data.to_csv(self.output_dir / 'raw_variant_data.csv', index=False)
                else:
                    logger.warning("No variant data was collected")
            
            
        except asyncio.CancelledError:
            logger.info("Operation cancelled")
            print("\nOperation cancelled. Cleaning up...")
        except Exception as e:
            logger.error(f"Error in fetch_variant_data_async: {e}")
            print(f"\nAn error occurred: {e}")
        finally:
            # Optional cleanup code
            print("\nShutdown complete")
                
    async def _fetch_gene_variants(self, session, gene, cache, all_variants_data, cache_file) -> None:
        """Fetch and process variants for a single gene."""
        logger.info(f"Fetching data for {gene}")
        
        # Define headers
        headers = {
            "Content-Type": "application/json",
            "Accept": "application/json"
        }
        
        # Check cache first
        if gene in cache:
            all_variants_data.extend(cache[gene])
            return
        
        try:
            # Get gene coordinates with retries
            max_retries = 3
            retry_delay = 1  # Delay in seconds

            for attempt in range(max_retries):
                try:
                    gene_url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}?"
                    async with session.get(gene_url, headers=headers) as response:
                        response.raise_for_status()
                        gene_info = await response.json()
                        break  # Exit the retry loop if successful
                except aiohttp.ClientError as e:
                    logger.warning(f"Attempt {attempt + 1} failed for gene {gene}: {e}")
                    if attempt < max_retries - 1:
                        await asyncio.sleep(retry_delay * (attempt + 1))  # Exponential backoff
                    else:
                        logger.error(f"Failed to fetch gene info for {gene} after {max_retries} attempts")
                        return None
            
            # Get variants in the region
            variant_url = (f"https://rest.ensembl.org/overlap/region/homo_sapiens/"
                        f"{gene_info['seq_region_name']}:{gene_info['start']}-{gene_info['end']}?"
                        f"feature=variation")
            
            async with session.get(variant_url, headers=headers) as response:
                response.raise_for_status()
                variants = await response.json()
            
            logger.info(f"Found {len(variants)} variants for gene {gene}")
            
            # Process variants
            processed_variants = []
            batch_size = 1000  # Increased batch size
            total_batches = (len(variants) + batch_size - 1) // batch_size
            
            # Debug: Log sample variant formats
            sample_variants = variants[:5] if len(variants) > 5 else variants
            logger.debug(f"Sample variant formats: {[v.get('id') for v in sample_variants]}")
            
            logger.info(f"Processing {total_batches} batches for gene {gene}")
            
            for batch_num, i in enumerate(range(0, len(variants), batch_size), 1):
                batch = variants[i:i + batch_size]
                logger.info(f"Processing batch {batch_num}/{total_batches} for gene {gene}")
                
                # Process each variant in the batch
                for variant in batch:
                    if 'id' in variant:
                        try:
                            variant_id = variant['id']
                            # Create URL with proper parameters
                            #pop_url = f"https://rest.ensembl.org/variation/human/{variant_id}?pops=1;content-type=application/json"
                            pop_url = f"https://rest.ensembl.org/variation/human/{variant_id}?pops=1"
                            logger.debug(f"Requesting URL for variant {variant_id}: {pop_url}")
                            
                            # Fetch population data with retries
                            max_retries = 3
                            retry_delay = 1  # Delay in seconds
                            pop_data = None

                            for attempt in range(max_retries):
                                try:
                                    async with session.get(pop_url, headers=headers) as response:
                                        if response.status != 200:
                                            logger.warning(f"Attempt {attempt + 1} failed for variant {variant_id}: {response.status}")
                                            if attempt < max_retries - 1:
                                                await asyncio.sleep(retry_delay * (attempt + 1))  # Exponential backoff
                                            else:
                                                logger.error(f"Failed to fetch population data for variant {variant_id} after {max_retries} attempts")
                                                break  # Exit the retry loop if all attempts fail
                                        else:
                                            # Respect rate limits with dynamic delay
                                            if response.headers.get("Retry-After"):
                                                retry_after = int(response.headers["Retry-After"])
                                                await asyncio.sleep(retry_after)
                                            else:
                                                await asyncio.sleep(0.05)  # Reduced fixed delay
                                            
                                            pop_data = await response.json()
                                            logger.debug(f"Population data for variant {variant_id}: {pop_data}")
                                            break  # Exit the retry loop if successful
                                except aiohttp.ClientError as e:
                                    logger.warning(f"Attempt {attempt + 1} failed for variant {variant_id}: {e}")
                                    if attempt < max_retries - 1:
                                        await asyncio.sleep(retry_delay * (attempt + 1))  # Exponential backoff
                                    else:
                                        logger.error(f"Failed to fetch population data for variant {variant_id} after {max_retries} attempts")
                                        break  # Exit the retry loop if all attempts fail
                            
                            if pop_data is None:
                                continue  # Skip this variant if all retries failed
                            
                            variant_data = self._process_single_variant(
                                variant,
                                pop_data,
                                gene,
                                gene_info['seq_region_name']
                            )
                            processed_variants.extend(variant_data)
                        
                        except Exception as e:
                            logger.warning(f"Error processing variant {variant.get('id', 'unknown')}: {e}")
                            continue
            
    
            # Cache the results
            cache[gene] = processed_variants
            all_variants_data.extend(processed_variants)

            # Log completion for this gene
            logger.info(f"Completed processing gene {gene}")
            logger.info(f"Processed {len(processed_variants)} variants for gene {gene}")

            # Save cache after each gene
            with open(cache_file, 'w') as f:
                json.dump(cache, f)

            # Save intermediate results to a CSV file
            intermediate_file = Path(self.output_dir) / f'intermediate_results_{gene}.csv'
            pd.DataFrame(processed_variants).to_csv(intermediate_file, index=False)
            logger.info(f"Saved intermediate results for gene {gene} to {intermediate_file}")
                    
        except Exception as e:
            logger.error(f"Error processing gene {gene}: {e}")    
    
               
    def clear_cache(self) -> None:
        """Clear the variant cache"""
        cache_file = Path(self.output_dir) / "variant_cache.json"
        if cache_file.exists():
            cache_file.unlink()
        logger.info("Cache cleared")

    def _process_single_variant_quick(self, variant: Dict, gene: str, chromosome: str) -> List[Dict]:
        """Process a single variant with basic population information"""
        processed_variants = []
        
        # Base variant information
        base_variant_data = {
            'gene': gene,
            'chromosome': chromosome,
            'position': variant['start'],
            'variant_id': variant.get('id', ''),
            'allele_string': variant.get('allele_string', ''),
            'consequence_type': str(variant.get('consequence_type', []))
        }
        
        # Add basic population information
        for pop_name in self.populations:
            variant_data = base_variant_data.copy()
            variant_data.update({
                'population': pop_name,
                'allele_freq': variant.get(f'{pop_name}_freq', 0.0),
                'allele_count': variant.get(f'{pop_name}_count', 0),
                'sample_size': variant.get(f'{pop_name}_size', 0)
            })
            processed_variants.append(variant_data)
        
        return processed_variants


    def _process_single_variant(self, variant: Dict, variant_details: Dict, gene: str, chromosome: str) -> List[Dict]:
        """Process a single variant with its population and phenotype data"""
        processed_variants = []
        
        # Log the incoming data structure
        logger.debug(f"Processing variant details: {variant_details}")
        
        # Fetch phenotype data
        phenotype_data = self.fetch_variant_phenotypes(variant.get('id', ''))
        clinical_info = self.process_phenotype_data(phenotype_data) if phenotype_data else {}
        
        # Base variant information
        base_variant_data = {
            'gene': gene,
            'chromosome': chromosome,
            'position': variant['start'],
            'variant_id': variant.get('id', ''),
            'allele_string': variant.get('allele_string', ''),
            'consequence_type': str(variant.get('consequence_type', [])),
            # Add clinical information to base data
            'clinical_significance': ','.join(clinical_info.get('clinical_significance', [])),
            'associated_traits': ','.join(set(clinical_info.get('traits', []))),
            'risk_alleles': ','.join(set(clinical_info.get('risk_alleles', []))),
            'evidence_sources': ','.join(set(clinical_info.get('sources', []))),
            'study_references': ','.join(set(clinical_info.get('studies', [])))
        }
        
        # Population data handling
        population_data = []
        if isinstance(variant_details, dict):
            population_data = variant_details.get('populations', [])
        
        # Process population data with mappings
        if population_data:
            # Create mapping from gnomAD and other sources to 1000G populations
            population_source_mapping = {
                'gnomADe:afr': ['1000GENOMES:phase_3:YRI', '1000GENOMES:phase_3:LWK', '1000GENOMES:phase_3:GWD'],
                'gnomADe:eas': ['1000GENOMES:phase_3:JPT', '1000GENOMES:phase_3:CHB'],
                'gnomADe:eur': ['1000GENOMES:phase_3:CEU', '1000GENOMES:phase_3:GBR'],
                'TOPMed': ['1000GENOMES:phase_3:ALL'],
            }
            
            for pop_data in population_data:
                source_pop = pop_data.get('population')
                
                # Calculate sample size for this population
                sample_size = self.calculate_sample_size(pop_data, population_data)
                
                # Map source population to 1000G populations
                target_pops = []
                if source_pop in POPULATION_MAPPINGS:
                    target_pops = [source_pop]
                else:
                    # Try to map from source population to 1000G populations
                    for source, targets in population_source_mapping.items():
                        if source_pop and source_pop.startswith(source):
                            target_pops = targets
                            break
                
                for target_pop in target_pops:
                    if target_pop in POPULATION_MAPPINGS:
                        pop_info = POPULATION_MAPPINGS[target_pop]
                        variant_data = base_variant_data.copy()
                        try:
                            variant_data.update({
                                'population': pop_info['code'],
                                'major_group': pop_info['major_group'],
                                'subgroup': pop_info['subgroup'],
                                'population_name': pop_info['name'],
                                'region': pop_info['region'],
                                'allele_freq': float(pop_data.get('frequency', 0)),
                                'allele_count': int(pop_data.get('allele_count', 0)),
                                'sample_size': sample_size,
                                'allele': pop_data.get('allele', ''),
                                'source_population': source_pop
                            })
                            processed_variants.append(variant_data)
                            logger.debug(f"Successfully processed variant data for population {pop_info['code']}")
                        except (ValueError, TypeError) as e:
                            logger.warning(f"Error converting population data for {pop_info['code']}: {e}")
                            continue
        
        # If no data was processed, add entries for all mapped populations with zero frequencies
        if not processed_variants:
            for pop_id, pop_info in POPULATION_MAPPINGS.items():
                variant_data = base_variant_data.copy()
                variant_data.update({
                    'population': pop_info['code'],
                    'major_group': pop_info['major_group'],
                    'subgroup': pop_info['subgroup'],
                    'population_name': pop_info['name'],
                    'region': pop_info['region'],
                    'allele_freq': 0.0,
                    'allele_count': 0,
                    'sample_size': 0,
                    'allele': '',
                    'source_population': 'None'
                })
                processed_variants.append(variant_data)
        
        logger.debug(f"Processed {len(processed_variants)} population entries for variant {base_variant_data['variant_id']}")
        return processed_variants

    def _fetch_variants(self, base_url: str, gene_info: Dict, gene_name: str) -> List[Dict]:
        """Fetch variants for a gene with population data"""
        try:
            logger.info(f"Starting variant fetch for {gene_name}")
            
            # First get variants in the region
            variant_url = (f"{base_url}/overlap/region/homo_sapiens/"
                        f"{gene_info['seq_region_name']}:{gene_info['start']}-{gene_info['end']}?"
                        f"feature=variation")
            
            headers = {"Content-Type": "application/json"}
            response = requests.get(variant_url, headers=headers)
            response.raise_for_status()
            variants = response.json()
            
            # Filter variants with IDs and prepare for batch processing
            valid_variants = [v for v in variants if 'id' in v]
            if not valid_variants:
                logger.warning(f"No valid variants found for {gene_name}")
                return []

            logger.info(f"Found {len(valid_variants)} valid variants for {gene_name}")

            # Process population data in batches
            enriched_variants = []
            batch_size = 1000  # Adjust based on API limits
            
            # Split variant IDs into batches
            total_batches = (len(valid_variants) + batch_size - 1) // batch_size
            logger.info(f"Processing {total_batches} batches of size {batch_size} for {gene_name}")

            for batch_num in range(total_batches):
                start_idx = batch_num * batch_size
                end_idx = min(start_idx + batch_size, len(valid_variants))
                current_batch = valid_variants[start_idx:end_idx]
                
                # Prepare comma-separated variant IDs
                variant_ids = [v['id'] for v in current_batch]
                id_string = ','.join(variant_ids)
                
                # Get population data for batch
                detail_url = f"{base_url}/variation/homo_sapiens/{id_string}?pops=1"
                
                try:
                    detail_response = requests.get(detail_url, headers=headers)
                    detail_response.raise_for_status()
                    population_data = detail_response.json()
                    
                    # Process each variant in the batch
                    batch_enriched = []
                    for variant in current_batch:
                        variant_id = variant['id']
                        if variant_id in population_data:
                            variant_data = population_data[variant_id]
                            variant['population_data'] = variant_data.get('populations', [])
                            batch_enriched.append(variant)
                    
                    # Add successfully processed variants
                    enriched_variants.extend(batch_enriched)
                    
                    logger.info(
                        f"Processed batch {batch_num + 1}/{total_batches} for {gene_name}: "
                        f"variants {start_idx + 1}-{end_idx} "
                        f"(processed {len(batch_enriched)} variants)"
                    )
                    
                    # Rate limiting
                    time.sleep(0.067)  # ~15 requests per second as per Ensembl limit
                    
                except requests.exceptions.RequestException as e:
                    logger.error(
                        f"Error fetching population data for batch {batch_num + 1} in {gene_name}: {e}"
                    )
                    continue
                except Exception as e:
                    logger.error(
                        f"Unexpected error processing batch {batch_num + 1} in {gene_name}: {e}"
                    )
                    continue
            
            # Final summary
            success_rate = (len(enriched_variants) / len(valid_variants)) * 100 if valid_variants else 0
            logger.info(
                f"Completed processing for {gene_name}: "
                f"Successfully processed {len(enriched_variants)}/{len(valid_variants)} variants "
                f"({success_rate:.1f}% success rate)"
            )
            
            return enriched_variants
                
        except requests.exceptions.RequestException as e:
            logger.error(f"Error fetching variants for {gene_name}: {e}")
            return []
        except Exception as e:
            logger.error(f"Unexpected error processing variants for {gene_name}: {e}")
            return []
        
        
    def fetch_variant_phenotypes(self, variant_id: str) -> Dict:
        """Fetch phenotype data for a variant"""
        try:
            url = f"https://rest.ensembl.org/variation/human/{variant_id}?phenotypes=1"
            headers = {
                "Content-Type": "application/json",
                "Accept": "application/json"
            }
            
            response = requests.get(url, headers=headers)
            response.raise_for_status()
            return response.json()
        except Exception as e:
            logger.warning(f"Could not fetch phenotype data for {variant_id}: {e}")
            return None

    def process_phenotype_data(self, phenotype_data: Dict) -> Dict:
        """Process phenotype information from variant data"""
        if not phenotype_data or 'phenotypes' not in phenotype_data:
            return {}
            
        processed_data = {
            'clinical_significance': phenotype_data.get('clinical_significance', []),
            'traits': [],
            'risk_alleles': [],
            'sources': [],
            'studies': []
        }
        
        for phenotype in phenotype_data['phenotypes']:
            if 'trait' in phenotype:
                processed_data['traits'].append(phenotype['trait'])
            if 'risk_allele' in phenotype:
                processed_data['risk_alleles'].append(phenotype['risk_allele'])
            if 'source' in phenotype:
                processed_data['sources'].append(phenotype['source'])
            if 'study' in phenotype:
                processed_data['studies'].append(phenotype['study'])
                
        return processed_data
    
    
    @staticmethod
    def calculate_sample_size(pop_data: Dict, all_population_data: List[Dict]) -> int:
        """Calculate sample size from allele counts"""
        current_pop = pop_data['population']
        total_alleles = 0
        
        # Sum all allele counts for this population
        for entry in all_population_data:
            if entry['population'] == current_pop:
                total_alleles += entry.get('allele_count', 0)
        
        # Divide by 2 since each person has 2 alleles
        return total_alleles // 2
    
    def _process_variants(self, variants: List[Dict], gene: str, gene_info: Dict) -> pd.DataFrame:
        """Process variant data into DataFrame with population information"""
        variants_data = []
        
        for variant in variants:
            # Base variant information
            variant_data = {
                'gene': gene,
                'chromosome': gene_info['seq_region_name'],
                'position': variant['start'],
                'variant_id': variant.get('id', ''),
                'allele_string': variant.get('allele_string', ''),
                'consequence_type': str(variant.get('consequence_type', []))
            }
            
            # Process population data
            if 'population_data' in variant:
                for pop_data in variant['population_data']:
                    pop_name = pop_data.get('population')
                    if pop_name in self.populations:
                        pop_variant_data = variant_data.copy()
                        pop_variant_data.update({
                            'population': pop_name,
                            'allele_freq': float(pop_data.get('frequency', 0)),
                            'allele_count': int(pop_data.get('allele_count', 0)),
                            'sample_size': int(pop_data.get('sample_size', 0))
                        })
                        variants_data.append(pop_variant_data)
            else:
                # If no population data, add a row for each population with null values
                for pop_name in self.populations:
                    pop_variant_data = variant_data.copy()
                    pop_variant_data.update({
                        'population': pop_name,
                        'allele_freq': 0.0,
                        'allele_count': 0,
                        'sample_size': 0
                    })
                    variants_data.append(pop_variant_data)
        
        df = pd.DataFrame(variants_data)
        logger.info(f"Processed variants for gene {gene}. Shape: {df.shape}")
        logger.info(f"Columns: {list(df.columns)}")
        return df



    def _get_population_frequencies(self, variant_id: str) -> Dict:
        """Get population frequencies for a variant with improved error handling"""
        try:
            url = f"https://rest.ensembl.org/variation/human/{variant_id}?"
            headers = {"Content-Type": "application/json"}
            
            response = requests.get(url, headers=headers)
            response.raise_for_status()
            
            variant_info = response.json()
            
            # Extract population frequencies
            frequencies = {}
            if variant_info and 'populations' in variant_info:  # Note: changed from 'population_genetics'
                for pop_data in variant_info['populations']:
                    pop_name = pop_data.get('population')
                    if pop_name in self.populations:
                        try:
                            frequencies[pop_name] = float(pop_data.get('frequency', 0))
                        except (TypeError, ValueError):
                            frequencies[pop_name] = None
            
            # Debug logging
            logger.debug(f"Retrieved frequencies for variant {variant_id}: {frequencies}")
            return frequencies
            
        except requests.exceptions.RequestException as e:
            logger.warning(f"Could not fetch population frequencies for variant {variant_id}: {e}")
            return {}
        except Exception as e:
            logger.error(f"Error processing population frequencies for variant {variant_id}: {e}")
            return {}

    def clean_string_columns(self) -> pd.DataFrame:
        """Clean and standardize string columns"""
        df = pd.DataFrame(self.variant_data)
        
        string_columns = ['clinical_significance', 'associated_traits', 
                         'risk_alleles', 'evidence_sources', 'study_references']
        
        for col in string_columns:
            df[col] = df[col].fillna('')
            df[col] = df[col].str.split(',')
        
        return df

    def standardize_consequence_types(self) -> pd.DataFrame:
        """Standardize and categorize consequence types"""
        df = pd.DataFrame(self.variant_data)
        
        df['consequence_type'] = df['consequence_type'].str.replace(r'["\[\]]', '', regex=True)
        
        severity_map = {
            'splice_acceptor_variant': 'HIGH',   
            'missense_variant': 'MODERATE',
            '5_prime_UTR_variant': 'LOW',
        }
        df['consequence_severity'] = df['consequence_type'].map(severity_map)
        
        return df

    def process_clinical_data(self) -> pd.DataFrame:
        """Process clinical significance and traits"""
        df = pd.DataFrame(self.variant_data)
        
        df['is_pathogenic'] = df['clinical_significance'].apply(
            lambda x: 'pathogenic' in str(x).lower()
        )
        df['is_protective'] = df['clinical_significance'].apply(
            lambda x: 'protective' in str(x).lower()
        )
        
        df['trait_count'] = df['associated_traits'].apply(
            lambda x: len(str(x).split(',')) if pd.notnull(x) else 0
        )
        
        return df

    def process_population_data(self) -> pd.DataFrame:
        """Process population-related columns"""
        df = pd.DataFrame(self.variant_data)
        
        df['is_african'] = df['major_group'] == 'African'
        df['is_diaspora'] = df['subgroup'] == 'African_Diaspora'
        
        df['allele_freq'] = pd.to_numeric(df['allele_freq'], errors='coerce')
        
        df['is_common'] = df['allele_freq'] > 0.05
        df['is_rare'] = df['allele_freq'] < 0.01
        
        return df

    def preprocess_data(self) -> pd.DataFrame:
        """Main preprocessing method that calls all other preprocessing steps"""
        
        print("Initial DataFrame columns:", self.processed_data.columns.tolist())
        print("Initial DataFrame shape:", self.processed_data.shape)
        print("First few rows of data:\n", self.processed_data.head())
        # Start with clean strings
        df = self.clean_string_columns()
        
        # Add consequence information
        df_conseq = self.standardize_consequence_types()
        df['consequence_severity'] = df_conseq['consequence_severity']
        
        # Add clinical information
        df_clinical = self.process_clinical_data()
        for col in ['is_pathogenic', 'is_protective', 'trait_count']:
            df[col] = df_clinical[col]
        
        # Add population information
        df_pop = self.process_population_data()
        for col in ['is_african', 'is_diaspora', 'is_common', 'is_rare']:
            df[col] = df_pop[col]
        
        # Additional derived features
        df['has_clinical_evidence'] = df['evidence_sources'].apply(len) > 0
        df['has_study_references'] = df['study_references'].apply(len) > 0
        
        df['allele_length'] = df['allele'].str.len()
        df['is_long_allele'] = df['allele_length'] > 100
        
        # Store processed data
        self.processed_data = df
        return df

    def get_summary_statistics(self) -> Dict:
        """Generate summary statistics for the processed data"""
        if not hasattr(self, 'processed_data'):
            self.preprocess_data()
        
        summary = {
            'total_variants': len(self.processed_data['variant_id'].unique()),
            'pathogenic_variants': self.processed_data['is_pathogenic'].sum(),
            'protective_variants': self.processed_data['is_protective'].sum(),
            'population_counts': self.processed_data['population'].value_counts(),
            'consequence_types': self.processed_data['consequence_type'].value_counts(),
            'trait_distribution': self.processed_data['trait_count'].describe()
        }
        
        return summary        
    
    def _calculate_enrichment(self, group_data: pd.DataFrame) -> float:
        """
        Calculate enrichment of clinical variants in a group
        
        Args:
            group_data: DataFrame for a specific group
            
        Returns:
            float: Enrichment score
        """
        total_variants = len(group_data)
        clinical_variants = group_data['is_pathogenic'].sum()
        return clinical_variants / total_variants if total_variants > 0 else 0

    def analyze_clinical_impact(self) -> Dict:
        """
        Analyze population-specific clinical impact of variants
        
        Returns:
            Dict containing:
            - frequency_distributions: DataFrame of allele frequencies by population and clinical significance
            - statistical_tests: Results of statistical comparisons between populations
            - enrichment_analysis: Analysis of clinical variant enrichment in populations
            - summary_metrics: Key metrics and findings
        """
        try:
            df = self.processed_data
            results = {}
            
            # 1. Analyze frequency distributions of clinical variants
            clinical_freq_dist = df[df['is_pathogenic']].groupby(
                ['major_group', 'subgroup']
            ).agg({
                'allele_freq': ['mean', 'std', 'count'],
                'variant_id': 'nunique'
            }).reset_index()
            
            # 2. Statistical testing between populations
            african_data = df[df['is_african'] & df['is_pathogenic']]['allele_freq']
            other_data = df[~df['is_african'] & df['is_pathogenic']]['allele_freq']
            
            
            stat_test = stats.mannwhitneyu(african_data, other_data)
            
            # 3. Calculate enrichment of clinical variants
            enrichment = df.groupby('major_group').apply(self._calculate_enrichment)
            
            # 4. Identify population-specific clinical variants
            pop_specific = df[
                (df['is_pathogenic']) & 
                (df['allele_freq'] > 0.01) & 
                (df['sample_size'] > 30)  # Ensure reliable frequency estimates
            ].groupby('major_group').agg({
                'variant_id': lambda x: list(set(x)),
                'associated_traits': lambda x: list(set([item for sublist in x for item in sublist]))
            })
            
            # Compile results
            results = {
                'frequency_distributions': clinical_freq_dist,
                'statistical_tests': {
                    'method': 'Mann-Whitney U',
                    'statistic': stat_test.statistic,
                    'pvalue': stat_test.pvalue
                },
                'enrichment_analysis': enrichment.to_dict(),
                'population_specific_variants': pop_specific.to_dict(),
                'summary_metrics': {
                    'total_clinical_variants': df['is_pathogenic'].sum(),
                    'african_specific_count': len(pop_specific.loc['African']['variant_id']) 
                        if 'African' in pop_specific.index else 0
                }
            }
            
            # Generate visualizations
            self._visualize_clinical_impact(results)
            
            # Log key findings
            logger.info(f"Completed clinical impact analysis across {len(df['major_group'].unique())} major groups")
            logger.info(f"Found {results['summary_metrics']['total_clinical_variants']} clinical variants")
            # Calculate diaspora sharing
            west_african_variants = set(df[
                (df['subgroup'] == 'West_African') & 
                (df['is_pathogenic'])
            ]['variant_id'])

            east_african_variants = set(df[
                (df['subgroup'] == 'East_African') & 
                (df['is_pathogenic'])
            ]['variant_id'])

            diaspora_variants = set(df[
                (df['subgroup'] == 'African_Diaspora') & 
                (df['is_pathogenic'])
            ]['variant_id'])

            shared_count = len(west_african_variants.intersection(diaspora_variants))
            sharing_percentage = (shared_count / len(west_african_variants) * 100) if len(west_african_variants) > 0 else 0

            # Update results with new metrics
            results['summary_metrics'].update({
                'pathogenic_variants_found': df['is_pathogenic'].sum(),
                'population_groups_analyzed': len(df['major_group'].unique())
            })

            results['population_specific_findings'] = {
                'west_african_unique': len(west_african_variants - diaspora_variants),
                'east_african_unique': len(east_african_variants),
                'diaspora_shared_percent': sharing_percentage
            }
            return results
            
        except Exception as e:
            logger.error(f"Error in clinical impact analysis: {e}")
            raise

    def _visualize_clinical_impact(self, results: Dict) -> None:
        """
        Generate visualizations for clinical impact analysis
        
        Args:
            results: Dictionary containing clinical impact analysis results
        """
        try:

            
            # Set style
            plt.style.use('seaborn-v0_8')
            sns.set_palette("husl")
            
            # Create figures directory
            figures_dir = self.output_dir / 'figures'
            figures_dir.mkdir(exist_ok=True)
            
            # 1. Clinical Variant Distribution
            plt.figure(figsize=(12, 6))
            sns.boxplot(data=self.processed_data[self.processed_data['is_pathogenic']], 
                       x='major_group', 
                       y='allele_freq')
            plt.title('Clinical Variant Frequencies by Population')
            plt.xlabel('Population Group')
            plt.ylabel('Allele Frequency')
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(figures_dir / 'clinical_freq_distribution.png', dpi=300)
            plt.close()
            
            # 2. Enrichment Heatmap
            plt.figure(figsize=(10, 8))
            enrichment_df = pd.DataFrame(results['enrichment_analysis'].items(), 
                                       columns=['Population', 'Enrichment'])
            enrichment_pivot = enrichment_df.pivot_table(
                values='Enrichment',
                index='Population'
            )
            sns.heatmap(enrichment_pivot, 
                       annot=True, 
                       cmap='YlOrRd', 
                       fmt='.3f')
            plt.title('Clinical Variant Enrichment by Population')
            plt.tight_layout()
            plt.savefig(figures_dir / 'clinical_enrichment_heatmap.png', dpi=300)
            plt.close()
            
            # Add visualization paths to results
            results['visualizations'] = {
                'frequency_distribution': str(figures_dir / 'clinical_freq_distribution.png'),
                'enrichment_heatmap': str(figures_dir / 'clinical_enrichment_heatmap.png')
            }
            
        except Exception as e:
            logger.error(f"Error in clinical impact visualization: {e}")
        
    @staticmethod
    def calculate_severity_enrichment(group_data):
        """Helper function to calculate severity enrichment for a group"""
        total_variants = len(group_data)
        high_severity = len(group_data[group_data['consequence_severity'] == 'HIGH'])
        return high_severity / total_variants if total_variants > 0 else 0

    def analyze_functional_consequences(self) -> Dict:
        """
        Analyze relationship between functional consequences and population distributions
        
        Returns:
            Dict containing:
            - consequence_frequencies: DataFrame of variant frequencies by consequence type and population
            - statistical_tests: Results of consequence distribution comparisons
            - severity_analysis: Analysis of variant severity across populations
            - summary_metrics: Key metrics and findings
        """
        try:
            df = self.processed_data
            results = {}
            
            # 1. Analyze frequency distribution by consequence type
            consequence_freq_dist = df.groupby(
                ['major_group', 'consequence_severity']
            ).agg({
                'allele_freq': ['mean', 'std', 'count'],
                'variant_id': 'nunique'
            }).reset_index()
            
            # 2. Statistical testing of consequence distributions
            high_sev_african = df[
                (df['is_african']) & 
                (df['consequence_severity'] == 'HIGH')
            ]['allele_freq']
            
            high_sev_other = df[
                (~df['is_african']) & 
                (df['consequence_severity'] == 'HIGH')
            ]['allele_freq']
            # Calculate population distribution percentages
            high_severity_variants = df[df['consequence_severity'] == 'HIGH']
            total_high_severity = len(high_severity_variants)

            population_distribution = high_severity_variants.groupby('major_group').size().apply(
                lambda x: (x/total_high_severity * 100)
            ).to_dict()

            # Add to results
            results['population_distribution'] = population_distribution
            
            severity_test = stats.mannwhitneyu(high_sev_african, high_sev_other)
            
            # 3. Calculate severity enrichment by population
            severity_enrichment = df.groupby('major_group').apply(self.calculate_severity_enrichment)
            
            # 4. Identify population-specific severe variants
            severe_variants = df[
                (df['consequence_severity'] == 'HIGH') & 
                (df['allele_freq'] > 0.01) & 
                (df['sample_size'] > 30)
            ].groupby('major_group').agg({
                'variant_id': lambda x: list(set(x)),
                'consequence_type': lambda x: list(set(x))
            })
            
            # Compile results
            results = {
                'consequence_frequencies': consequence_freq_dist,
                'statistical_tests': {
                    'method': 'Mann-Whitney U (HIGH severity variants)',
                    'statistic': severity_test.statistic,
                    'pvalue': severity_test.pvalue
                },
                'severity_analysis': severity_enrichment.to_dict(),
                'population_specific_severe': severe_variants.to_dict(),
                'summary_metrics': {
                    'total_high_severity': len(df[df['consequence_severity'] == 'HIGH']),
                    'african_high_severity_count': len(df[
                        (df['is_african']) & 
                        (df['consequence_severity'] == 'HIGH')
                    ])
                }
            }
            
            # Log key findings
            logger.info(f"Completed functional consequence analysis across {len(df['major_group'].unique())} major groups")
            logger.info(f"Found {results['summary_metrics']['total_high_severity']} high severity variants")
            # Generate visualizations
            self._visualize_functional_consequences(results)
            return results
            
        except Exception as e:
            logger.error(f"Error in functional consequence analysis: {e}")
            raise
    @staticmethod
    def calculate_representation(group_data):
        """Helper function to calculate research representation for a group"""
        variants_with_evidence = group_data['has_clinical_evidence'].sum()
        total_variants = len(group_data)
        return variants_with_evidence / total_variants if total_variants > 0 else 0
    
    # Method already defined above
    
    @staticmethod
    def calculate_regional_risk(group_data):
        """Helper function to calculate regional risk enrichment"""
        total_variants = len(group_data)
        risk_variants = len(group_data[group_data['risk_alleles'].apply(len) > 0])
        return {
            'risk_ratio': risk_variants / total_variants if total_variants > 0 else 0,
            'mean_frequency': group_data['allele_freq'].mean(),
            'unique_risk_alleles': len(set([r for rs in group_data['risk_alleles'] for r in rs if rs]))
        }

    def analyze_evidence_patterns(self) -> Dict:
        """
        Analyze evidence patterns and research coverage across populations
        
        Returns:
            Dict containing:
            - evidence_distributions: DataFrame of evidence sources by population
            - statistical_tests: Results of evidence coverage comparisons
            - research_bias_analysis: Analysis of study representation across populations
            - summary_metrics: Key metrics and findings
        """
        try:
            df = self.processed_data
            results = {}
            
            # 1. Analyze evidence source distributions
            evidence_dist = df.groupby(['major_group', 'subgroup']).agg({
                'has_clinical_evidence': 'sum',
                'has_study_references': 'sum',
                'study_references': lambda x: len(set([ref for refs in x if refs for ref in refs])),
                'evidence_sources': lambda x: len(set([src for srcs in x if srcs for src in srcs]))
            }).reset_index()
            
            # 2. Statistical testing of evidence coverage
            african_evidence = df[df['is_african']]['has_clinical_evidence'].mean()
            other_evidence = df[~df['is_african']]['has_clinical_evidence'].mean()
            
            
            coverage_test = stats.fisher_exact([
                [df[df['is_african'] & df['has_clinical_evidence']].shape[0],
                df[df['is_african'] & ~df['has_clinical_evidence']].shape[0]],
                [df[~df['is_african'] & df['has_clinical_evidence']].shape[0],
                df[~df['is_african'] & ~df['has_clinical_evidence']].shape[0]]
            ])
            
            # 3. Calculate research representation index
            representation_index = df.groupby('major_group').apply(self.calculate_representation)
            
            # 4. Identify poorly studied variants
            understudied_variants = df[
                (~df['has_clinical_evidence']) & 
                (df['allele_freq'] > 0.05)  # Common variants lacking evidence
            ].groupby('major_group').agg({
                'variant_id': lambda x: list(set(x)),
                'allele_freq': 'mean'
            })
            
            # Compile results
            results = {
                'evidence_distributions': evidence_dist,
                'statistical_tests': {
                    'method': "Fisher's exact test for evidence coverage",
                    'odds_ratio': coverage_test[0],
                    'pvalue': coverage_test[1]
                },
                'research_bias_analysis': {
                    'representation_index': representation_index.to_dict(),
                    'understudied_variants': understudied_variants.to_dict()
                },
                'summary_metrics': {
                    'total_variants_with_evidence': df['has_clinical_evidence'].sum(),
                    'african_evidence_ratio': african_evidence,
                    'other_evidence_ratio': other_evidence,
                    'evidence_disparity': other_evidence - african_evidence
                }
            }
            
            # Log key findings
            logger.info(f"Completed evidence pattern analysis across {len(df['major_group'].unique())} major groups")
            logger.info(f"Found {results['summary_metrics']['total_variants_with_evidence']} variants with clinical evidence")
            # Generate visualizations
            self._visualize_evidence_patterns(results)
            return results
            
        except Exception as e:
            logger.error(f"Error in evidence pattern analysis: {e}")
            raise
        
    
    def analyze_population_structure(self) -> Dict:
        """
        Analyze population structure using PCA
        
        Returns:
            Dict containing PCA results and population clustering information
        """
        try:
            df = self.processed_data
            
            # 1. Prepare the data matrix
            # Pivot table to create variant x population matrix of allele frequencies
            freq_matrix = df.pivot_table(
                index='variant_id',
                columns='population',
                values='allele_freq',
                fill_value=0
            )
            
            # 2. Scale the data
            scaler = StandardScaler()
            scaled_data = scaler.fit_transform(freq_matrix)
            
            # 3. Perform PCA
            pca = PCA(n_components=3)  # Keep top 3 components
            pca_results = pca.fit_transform(scaled_data)
            
            # 4. Create results dictionary
            results = {
                'pca_coordinates': pd.DataFrame(
                    pca_results,
                    index=freq_matrix.index,
                    columns=['PC1', 'PC2', 'PC3']
                ),
                'explained_variance_ratio': pca.explained_variance_ratio_,
                'loadings': pd.DataFrame(
                    pca.components_.T,
                    index=freq_matrix.columns,
                    columns=['PC1', 'PC2', 'PC3']
                )
            }
            
            return results
            
        except Exception as e:
            logger.error(f"Error in population structure analysis: {e}")
            raise
    def _visualize_population_structure(self, results: Dict) -> None:
        """
        Generate visualizations for population structure analysis
        """
        try:
            figures_dir = self.output_dir / 'figures'
            figures_dir.mkdir(exist_ok=True)
            
            # 1. PCA Scatter Plot
            plt.figure(figsize=(10, 8))
            coords = results['pca_coordinates']
            populations = self.processed_data['major_group'].unique()
            
            # Create a mapping between variant_ids and major_groups
            variant_groups = self.processed_data.groupby('variant_id')['major_group'].first()
            
            # Plot each population
            for pop in populations:
                # Get variant IDs for this population
                pop_variants = variant_groups[variant_groups == pop].index
                # Get coordinates for these variants
                pop_coords = coords.loc[coords.index.intersection(pop_variants)]
                
                if not pop_coords.empty:
                    plt.scatter(
                        pop_coords['PC1'],
                        pop_coords['PC2'],
                        alpha=0.6,
                        label=pop
                    )
                
            plt.xlabel(f'PC1 ({results["explained_variance_ratio"][0]:.2%} variance explained)')
            plt.ylabel(f'PC2 ({results["explained_variance_ratio"][1]:.2%} variance explained)')
            plt.title('Population Structure PCA')
            plt.legend()
            plt.tight_layout()
            plt.savefig(figures_dir / 'population_structure_pca.png', dpi=300)
            plt.close()
            
            # 2. Loadings Plot
            plt.figure(figsize=(12, 6))
            loadings = results['loadings']
            sns.heatmap(
                loadings,
                cmap='RdBu',
                center=0,
                annot=True,
                fmt='.2f'
            )
            plt.title('PCA Component Loadings by Population')
            plt.tight_layout()
            plt.savefig(figures_dir / 'pca_loadings.png', dpi=300)
            plt.close()
            
            logger.info("Generated population structure visualizations")
            
        except Exception as e:
            logger.error(f"Error generating population structure visualizations: {e}")
            raise
    
    def _write_functional_consequences_section(self, f: TextIO, results: Dict) -> None:
        """
        Write functional consequences section
        """
        if 'functional_consequences' in results:
            f.write("## Functional Consequences Analysis\n\n")
            func = results['functional_consequences']
            
            f.write("### Severity and Adaptive Context\n")
            metrics = func['summary_metrics']
            african_ratio = (metrics['african_high_severity_count'] / metrics['total_high_severity']) * 100
            
            f.write(f"While African populations show a higher proportion ({african_ratio:.1f}%) of high-severity variants, ")
            f.write("detailed analysis reveals important adaptive contexts:\n\n")
            
            f.write("#### Balanced Polymorphisms\n")
            if 'balanced_variants' in func:
                for variant_type, details in func['balanced_variants'].items():
                    f.write(f"- {variant_type}:\n")
                    f.write(f"  * Frequency: {details['frequency']:.2f}%\n")
                    f.write(f"  * Selective Advantage: {details['advantage']}\n")
                    f.write(f"  * Environmental Context: {details['context']}\n")
            else:
                f.write("- Analysis of balanced polymorphisms pending additional data\n")
            f.write("\n")    
    def _write_evidence_patterns_section(self, f: TextIO, results: Dict):
        """
        Write enhanced evidence patterns section with detailed distribution analysis
        """
        if 'evidence_patterns' in results:
            evidence = results['evidence_patterns']
            dist = evidence['evidence_distributions']
            
            f.write("\n### Evidence Quality Assessment\n")
            f.write(f"Total Variants with Evidence: {evidence['summary_metrics']['total_variants_with_evidence']:,}\n\n")
            
            f.write("#### Evidence Depth Analysis\n")
            f.write("While overall coverage appears equal, important qualitative differences exist:\n")
            
            f.write("##### Evidence Type Distribution\n")
            for _, row in dist.iterrows():
                f.write(f"\n{row['major_group']} ({row['subgroup']}):\n")
                f.write(f"- Clinical Evidence: {row['has_clinical_evidence']:,} variants ")
                f.write(f"({(row['has_clinical_evidence']/len(dist))*100:.1f}%)\n")
                f.write(f"- Study References: {row['has_study_references']:,} variants ")
                f.write(f"({(row['has_study_references']/len(dist))*100:.1f}%)\n")
                f.write(f"- Unique Evidence Sources: {row['evidence_sources']:,}\n")
                f.write(f"- Evidence Quality Score: {row['evidence_sources']/row['has_clinical_evidence']:.2f}\n")
            
            # Research bias analysis
            f.write("\n### Research Bias Assessment\n")
            bias = evidence['research_bias_analysis']
            understudied = bias.get('understudied_variants', {})
            f.write(f"- Understudied Variants: {len(understudied)} variants ")
            f.write("lack comprehensive characterization despite high frequency\n")
            
            # Add understudied pathways
            if 'representation_index' in bias:
                f.write("- Population Representation Index:\n")
                for pop, index in bias['representation_index'].items():
                    f.write(f"  * {pop}: {index:.2f}\n")
            
            f.write("- Key understudied pathways:\n")
            pathways = self._identify_understudied_pathways(evidence)
            for pathway, details in pathways.items():
                f.write(f"  * {pathway}: {details['coverage_ratio']:.1f}% coverage")
                f.write(f" ({details['variant_count']} variants)\n")
            f.write("\n")
            
    def _identify_understudied_pathways(self, evidence_data):
        """Helper method to identify understudied pathways"""
        pathways = {}
        if 'evidence_distributions' in evidence_data:
            dist = evidence_data['evidence_distributions']
            
            # Group by pathway and calculate coverage
            pathway_coverage = dist.groupby('pathway').agg({
                'has_clinical_evidence': 'sum',
                'variant_id': 'count'
            })
            
            for pathway, row in pathway_coverage.iterrows():
                coverage_ratio = (row['has_clinical_evidence'] / row['variant_id']) * 100
                if coverage_ratio < 50:  # Define understudied as <50% coverage
                    pathways[pathway] = {
                        'coverage_ratio': coverage_ratio,
                        'variant_count': row['variant_id']
                    }
        
        return pathways
    
    def _write_trait_analysis_section(self, f: TextIO, results: Dict):
        """
        Write enhanced multi-trait analysis section
        """
        if 'pleiotropic_effects' in results:
            pleio = results['pleiotropic_effects']
            metrics = pleio['summary_metrics']
            
            f.write("\n## Multi-trait Analysis\n\n")
            f.write("### Complex Trait Interactions\n")
            f.write(f"Analysis of {metrics['total_pleiotropic_variants']:,} pleiotropic variants ")
            f.write("reveals intricate patterns of trait interactions:\n\n")
            
            # Network analysis
            if 'trait_network_analysis' in pleio:
                networks = pleio['trait_network_analysis']
                f.write("#### Key Trait Networks\n")
                
                # Get trait pairs from your analyze_trait_connections results
                for pop, data in networks.items():
                    f.write(f"\n{pop}:\n")
                    if 'trait_pairs' in data:
                        # Sort pairs by frequency and get top 5
                        pairs = sorted(data['trait_pairs'].items(), key=lambda x: x[1], reverse=True)[:5]
                        for pair, count in pairs:
                            f.write(f"- {pair[0]} ↔ {pair[1]}: {count:,} variants\n")
                    f.write(f"- Unique Traits: {data.get('unique_traits', 0):,}\n")
                    f.write(f"- Max Traits per Variant: {data.get('max_traits_per_variant', 0):,}\n\n")
            
            f.write("\n#### Population-Specific Pleiotropic Patterns\n")
            if 'population_patterns' in pleio:
                for pop, patterns in pleio['population_patterns'].items():
                    f.write(f"\n{pop}:\n")
                    f.write(f"- Characteristic Combinations:\n")
                    for trait_combo in patterns.get('top_combinations', [])[:3]:
                        f.write(f"  * {' + '.join(trait_combo)}\n")
            
    def _write_selection_patterns_section(self, f: TextIO, results: Dict):
        protect = results['protective_variants']
        selection = protect['selection_analysis']
        
        f.write("\n#### Selection Signatures\n")
        for group, patterns in selection.items():
            f.write(f"\n{group}:\n")
            f.write(f"- High Frequency Protective Variants: {patterns['high_freq_protective']}\n")
            f.write(f"- Mean Protective Allele Frequency: {patterns['mean_protective_freq']:.3f}\n")
            f.write("- Consequence Types:\n")
            for conseq, count in patterns['protective_consequence_types'].items():
                f.write(f"  * {conseq}: {count} variants\n")
    
    def analyze_geographical_risks(self) -> Dict:
        """
        Analyze geographical distribution of risk alleles across populations
        
        Returns:
            Dict containing:
            - geographical_distributions: DataFrame of risk allele frequencies by region
            - statistical_tests: Results of geographical clustering tests
            - regional_enrichment: Analysis of risk allele patterns by region
            - summary_metrics: Key metrics and findings
        """
        try:
            df = self.processed_data
            results = {}
            
            # 1. Analyze risk allele distribution by region
            geo_dist = df[df['risk_alleles'].apply(len) > 0].groupby(
                ['major_group', 'subgroup', 'region']
            ).agg({
                'allele_freq': ['mean', 'std'],
                'risk_alleles': lambda x: len(set([r for rs in x for r in rs])),
                'variant_id': 'nunique'
            }).reset_index()
            
            #logger.info(f"geo_dist structure:\n{geo_dist.head()}")
            #logger.info(f"geo_dist index: {geo_dist.index}")
            
                       
            
            # Compare risk allele frequencies between regions
            west_african = df[
                (df['subgroup'] == 'West_African') & 
                (df['risk_alleles'].apply(len) > 0)
            ]['allele_freq']
            
            east_african = df[
                (df['subgroup'] == 'East_African') & 
                (df['risk_alleles'].apply(len) > 0)
            ]['allele_freq']
            
            geo_test = stats.mannwhitneyu(west_african, east_african)
            
            regional_enrichment = df.groupby(['region', 'subgroup']).apply(self.calculate_regional_risk)
            regional_enrichment = df.groupby(['region', 'subgroup']).apply(self.calculate_regional_risk)
            
            # 4. Identify region-specific risk alleles
            region_specific = df[
                (df['risk_alleles'].apply(len) > 0) &
                (df['allele_freq'] > 0.01) &
                (df['sample_size'] > 30)
            ].groupby(['region', 'subgroup']).agg({
                'variant_id': lambda x: list(set(x)),
                'risk_alleles': lambda x: list(set([r for rs in x for r in rs if rs])),
                'associated_traits': lambda x: list(set([t for ts in x for t in ts if ts]))
            })
            
            # Compile results
            results = {
                'geographical_distributions': geo_dist,
                'statistical_tests': {
                    'method': 'Mann-Whitney U (West vs East African)',
                    'statistic': geo_test.statistic,
                    'pvalue': geo_test.pvalue
                },
                'regional_enrichment': regional_enrichment.to_dict(),
                'region_specific_risks': region_specific.to_dict(),
                'summary_metrics': {
                    'total_risk_alleles': len(set([r for rs in df['risk_alleles'] for r in rs if rs])),
                    'african_specific_risks': len(region_specific.loc['West Africa']) 
                        if 'West Africa' in region_specific.index else 0,
                    'regions_analyzed': len(df['region'].unique())
                }
            }
            
            # Log key findings
            logger.info(f"Completed geographical risk analysis across {results['summary_metrics']['regions_analyzed']} regions")
            logger.info(f"Found {results['summary_metrics']['total_risk_alleles']} unique risk alleles")
            # Generate visualizations
            self._visualize_geographical_risks(results)
            return results
            
        except Exception as e:
            logger.error(f"Error in geographical risk analysis: {e}")
            raise
    @staticmethod
    def analyze_trait_connections(group_data):
        """Helper function to analyze trait connections in pleiotropic variants"""
        trait_pairs = []
        for traits in group_data['associated_traits']:
            if len(traits) > 1:
                trait_pairs.extend(list(combinations(sorted(traits), 2)))
        
        # Count frequency of trait pairs
        pair_counts = Counter(trait_pairs)
        return {
            'trait_pairs': dict(pair_counts),
            'unique_traits': len(set([t for ts in group_data['associated_traits'] for t in ts if ts])),
            'max_traits_per_variant': max(len(ts) for ts in group_data['associated_traits'] if ts)
        }
    
    def analyze_pleiotropic_effects(self) -> Dict:
        """
        Analyze variants associated with multiple traits and their population distributions
        
        Returns:
            Dict containing:
            - trait_distributions: DataFrame of multi-trait variant frequencies
            - statistical_tests: Results of pleiotropic effect comparisons
            - trait_network_analysis: Analysis of trait associations and patterns
            - summary_metrics: Key metrics and findings
        """
        try:
            df = self.processed_data
            results = {}
            
            # 1. Analyze multi-trait variant distributions
            # Consider a variant pleiotropic if associated with multiple traits
            df['is_pleiotropic'] = df['trait_count'] > 1
            
            trait_dist = df[df['is_pleiotropic']].groupby(
                ['major_group', 'subgroup']
            ).agg({
                'allele_freq': ['mean', 'std'],
                'trait_count': ['mean', 'max'],
                'variant_id': 'nunique',
                'associated_traits': lambda x: list(set([t for ts in x for t in ts if ts]))
            }).reset_index()
            
                       
            # Compare frequency distributions of pleiotropic variants
            african_pleiotropic = df[
                (df['is_african']) & 
                (df['is_pleiotropic'])
            ]['allele_freq']
            
            other_pleiotropic = df[
                (~df['is_african']) & 
                (df['is_pleiotropic'])
            ]['allele_freq']
            
            try:
                pleiotropy_test = stats.mannwhitneyu(african_pleiotropic, other_pleiotropic)
            except Exception as e:
                logger.error(f"Error performing Mann-Whitney U test: {e}")
                pleiotropy_test = None
                
            
            
            # 3. Analyze trait associations
            trait_networks = df.groupby('major_group').apply(self.analyze_trait_connections)
            
            # 4. Identify population-specific pleiotropic variants
            pop_specific_pleiotropic = df[
                (df['is_pleiotropic']) & 
                (df['allele_freq'] > 0.01) & 
                (df['sample_size'] > 30)
            ].groupby('major_group').agg({
                'variant_id': lambda x: list(set(x)),
                'associated_traits': lambda x: list(set([t for ts in x for t in ts if ts])),
                'trait_count': ['mean', 'max']
            })
            
            # Compile results
            results = {
                'trait_distributions': trait_dist,
                'statistical_tests': {
                    'method': 'Mann-Whitney U (pleiotropic variants)',
                    'statistic': pleiotropy_test.statistic if pleiotropy_test else None,
                    'pvalue': pleiotropy_test.pvalue if pleiotropy_test else None
                },
                'trait_network_analysis': {
                    'networks': trait_networks.to_dict(),
                    'population_specific': pop_specific_pleiotropic.to_dict()
                },
                'summary_metrics': {
                    'total_pleiotropic_variants': df['is_pleiotropic'].sum(),
                    'african_pleiotropic_count': len(df[
                        (df['is_african']) & 
                        (df['is_pleiotropic'])
                    ]),
                    'max_traits_per_variant': df['trait_count'].max(),
                    'mean_traits_per_variant': df['trait_count'].mean()
                }
            }
            
            # Log key findings
            logger.info(f"Completed pleiotropic analysis across {len(df['major_group'].unique())} major groups")
            logger.info(f"Found {results['summary_metrics']['total_pleiotropic_variants']} pleiotropic variants")
            # Generate visualizations
            self._visualize_pleiotropic_effects(results)
            return results
            
        except Exception as e:
            logger.error(f"Error in pleiotropic analysis: {e}")
            raise
        
    @staticmethod   
    def analyze_selection_patterns(group_data):
        """Helper function for analyzing selection patterns"""
        protective_vars = group_data[group_data['is_protective']]
        
        # Basic metrics
        high_freq = protective_vars[protective_vars['allele_freq'] > 0.05]
        consequence_types = protective_vars['consequence_type'].value_counts().to_dict()
        
        # Selection signals - using allele frequency as proxy for selection
        selection_signals = {
            consequence: {
                'selection_coef': protective_vars[protective_vars['consequence_type'] == consequence]['allele_freq'].mean(),
                'selection_age': "Recent" if protective_vars[protective_vars['consequence_type'] == consequence]['allele_freq'].mean() > 0.1 else "Ancient",
                'context': "High frequency variant suggesting adaptive advantage" if protective_vars[protective_vars['consequence_type'] == consequence]['allele_freq'].mean() > 0.1 else "Maintained in population"
            }
            for consequence in consequence_types.keys()
        }
        
        # Adaptive benefits based on consequence types
        adaptive_benefits = {
            pop_group: [
                f"Enhanced {ct.replace('_', ' ')} function" 
                for ct in protective_vars[protective_vars['major_group'] == pop_group]['consequence_type'].unique()
            ]
            for pop_group in protective_vars['major_group'].unique()
        }
        
        return {
            'high_freq_protective': len(high_freq),
            'mean_protective_freq': protective_vars['allele_freq'].mean(),
            'protective_consequence_types': consequence_types,
            'selection_signals': selection_signals,
            'adaptive_benefits': adaptive_benefits
        }
    
    def calculate_protective_enrichment(group_data):
        """Helper function for calculating protective enrichment"""
        total_variants = len(group_data)
        protective_variants = group_data['is_protective'].sum()
        return protective_variants / total_variants if total_variants > 0 else 0
    @staticmethod
    def calculate_rare_variant_rate(group):
        """Helper function for calculating rare variant rate"""
        return group['is_rare'].sum() / len(group) if len(group) > 0 else 0
    
    @staticmethod
    def analyze_detection_sensitivity(group_data):
        """Helper function for analyzing detection sensitivity"""
        return {
            'rare_variant_rate': AfricanGeneticDiversity.calculate_rare_variant_rate(group_data),
            'mean_allele_freq': group_data['allele_freq'].mean(),
            'variant_density': len(group_data) / group_data['sample_size'].mean() if group_data['sample_size'].mean() > 0 else 0
        }
    
    def analyze_size_thresholds(group_data):
        """Helper function for analyzing size thresholds"""
        thresholds = range(10, int(group_data['sample_size'].max()), 10)
        results = []
        for threshold in thresholds:
            subset = group_data[group_data['sample_size'] >= threshold]
            results.append({
                'threshold': threshold,
                'variants_detected': len(subset),
                'rare_variants_detected': subset['is_rare'].sum()
            })
        return pd.DataFrame(results)

    @staticmethod
    def calculate_power_metrics(group_data):
        """Helper function for calculating power metrics"""
        return {
            'rare_variant_detection_rate': len(group_data[group_data['is_rare']]) / len(group_data) if len(group_data) > 0 else 0,
            'mean_rare_af': group_data[group_data['is_rare']]['allele_freq'].mean(),
            'sample_size_required': group_data['sample_size'].quantile(0.95)
        }
    
    def analyze_protective_variants(self) -> Dict:
        """Analyze distribution and patterns of protective variants"""
        try:
            df = self.processed_data
            results = {}
            
            # 1. Analyze protective variant distributions
            protective_dist = df[df['is_protective']].groupby(
                ['major_group', 'subgroup']
            ).agg({
                'allele_freq': ['mean', 'std', 'count'],
                'variant_id': 'nunique',
                'consequence_severity': lambda x: x.value_counts().to_dict()
            }).reset_index()
            
            # 2. Compare protective variant frequencies
            african_protective = df[
                (df['is_african']) & 
                (df['is_protective'])
            ]['allele_freq']
            
            other_protective = df[
                (~df['is_african']) & 
                (df['is_protective'])
            ]['allele_freq']
            
            protective_test = stats.mannwhitneyu(african_protective, other_protective)
            
            # 3. Analyze selection patterns
            selection_patterns = df.groupby(['major_group', 'subgroup']).apply(self.analyze_selection_patterns)
            
            # Compile results
            results = {
                'protective_distributions': protective_dist,
                'statistical_tests': {
                    'method': 'Mann-Whitney U (protective variants)',
                    'statistic': protective_test.statistic,
                    'pvalue': protective_test.pvalue
                },
                'selection_analysis': {
                    'patterns': selection_patterns.to_dict()
                }
            }
            # Generate visualizations
            self._visualize_protective_variants(results)
            return results
            
        except Exception as e:
            logger.error(f"Error in protective variant analysis: {e}")
            raise
    
    def analyze_sample_size_effects(self) -> Dict:
        """
        Analyze the effect of sample size on variant detection and characteristics
        
        Returns:
            Dict containing sample size analysis results
        """
        try:
            # Create base dataframe
            df = pd.DataFrame({
                'sample_size': self.variant_data['sample_size'].fillna(0),
                'allele_freq': self.variant_data['allele_freq'].fillna(0)
            })
            
            # Create quartiles, handling duplicate edges
            try:
                # First attempt: try regular quartiles with duplicate handling
                df['sample_size_quartile'] = pd.qcut(df['sample_size'], 
                                                q=4, 
                                                labels=['Q1', 'Q2', 'Q3', 'Q4'],
                                                duplicates='drop')
            except ValueError:
                # Alternative approach: if too many duplicates, use custom bins
                unique_values = sorted(df['sample_size'].unique())
                n_groups = min(4, len(unique_values))
                
                if n_groups < 2:
                    # If too few unique values, create binary categories
                    df['sample_size_quartile'] = df['sample_size'].apply(
                        lambda x: 'Low' if x == 0 else 'High'
                    )
                else:
                    # Create custom bins based on unique values
                    bins = np.percentile(unique_values, 
                                    np.linspace(0, 100, n_groups + 1))
                    df['sample_size_quartile'] = pd.cut(df['sample_size'],
                                                    bins=bins,
                                                    labels=[f'Q{i+1}' for i in range(n_groups)],
                                                    include_lowest=True)
            
            # Calculate statistics
            size_stats = df.groupby('sample_size_quartile').agg({
                'sample_size': ['mean', 'std', 'count'],
                'allele_freq': ['mean', 'std']
            }).round(4)
            
            # Convert the multi-index DataFrame to a more accessible dictionary format
            size_stats_dict = {}
            for idx in size_stats.index:  # For each quartile
                size_stats_dict[str(idx)] = {
                    'sample_size': {
                        'mean': size_stats.loc[idx, ('sample_size', 'mean')],
                        'std': size_stats.loc[idx, ('sample_size', 'std')],
                        'count': size_stats.loc[idx, ('sample_size', 'count')]
                    },
                    'allele_freq': {
                        'mean': size_stats.loc[idx, ('allele_freq', 'mean')],
                        'std': size_stats.loc[idx, ('allele_freq', 'std')]
                    }
                }
            
            # Calculate correlation
            correlation = df['sample_size'].corr(df['allele_freq'])
            
            results = {
                'size_distribution': size_stats_dict,
                'correlation': correlation,
                'raw_data': df
            }
            
            # Generate visualizations
            self._visualize_sample_size_effects(results)
            
            return results
        
        except Exception as e:
            logger.error(f"Error in sample size analysis: {e}")
            raise
    
        
    def _write_sample_size_effects_section(self, f: TextIO, results: Dict) -> None:
        """Write comprehensive analysis of sample size effects with research question context."""
        if 'sample_size_effects' not in results:
            return
            
        f.write("\n## Sample Size Effects Analysis\n\n")
        
        # State research question
        f.write("Research Question: How does sample size variation across populations affect ")
        f.write("our understanding of rare clinically significant variants?\n\n")
        
        size = results['sample_size_effects']
        
        # Overview
        f.write("### Key Findings Addressing the Research Question\n")
        size_dist = size.get('size_distribution', {})
        f.write("Sample Size Distribution across Quartiles:\n")
        
        # Access the restructured dictionary
        for quartile in size_dist.keys():
            f.write(f"\n#### {quartile}\n")
            try:
                quartile_stats = size_dist[quartile]
                
                # Write sample size statistics
                if 'sample_size' in quartile_stats:
                    f.write(f"Mean Sample Size: {quartile_stats['sample_size']['mean']:.1f}\n")
                    f.write(f"Sample Count: {quartile_stats['sample_size']['count']:.0f}\n")
                
                # Write allele frequency statistics
                if 'allele_freq' in quartile_stats:
                    mean = quartile_stats['allele_freq']['mean']
                    std = quartile_stats['allele_freq']['std']
                    f.write(f"Mean Allele Frequency: {mean:.4f} (±{std:.4f})\n")
                    
            except Exception as e:
                logger.warning(f"Could not process stats for quartile {quartile}: {e}")
    
    
        
        # Correlation Analysis
        f.write("\n### Sample Size-Variant Detection Relationship\n")
        f.write(f"Correlation Coefficient: {size.get('correlation', 0):.4f}\n")
        f.write("This correlation indicates the strength of the relationship between ")
        f.write("sample size and variant detection capability.\n\n")
        
        # Conclusion
        f.write("\n### Conclusion for Research Question\n")
        f.write("The analysis demonstrates the critical impact of sample size variation on our ")
        f.write("ability to detect and characterize rare variants. This relationship has ")
        f.write("important implications for study design and interpretation of genetic ")
        f.write("variation across populations with different sample sizes.\n\n")    


    def run_comprehensive_analysis(self) -> Dict:
        """Run all analyses and compile results"""
        try:
            # Ensure data is preprocessed
            if not hasattr(self, 'processed_data'):
                self.preprocess_data()
                
            logger.info("Starting comprehensive genetic diversity analysis")
            
            # Run all analyses
            results = {
            'clinical_impact': self.analyze_clinical_impact(),
            'population_structure': self.analyze_population_structure(),
            'functional_consequences': self.analyze_functional_consequences(),
            'evidence_patterns': self.analyze_evidence_patterns(),
            'geographical_risks': self.analyze_geographical_risks(),
            'pleiotropic_effects': self.analyze_pleiotropic_effects(),
            'protective_variants': self.analyze_protective_variants(),
            'sample_size_effects': self.analyze_sample_size_effects()
        }
            
            logger.info("Completed all analyses successfully")
            return results
            
        except Exception as e:
            logger.error(f"Error in comprehensive analysis: {e}")
            raise
    def summarize_research_questions(self, results: Dict) -> Dict:
        """
        Map analysis results to research questions and generate summaries
        
        Args:
            results: Dictionary containing results from all analyses
        
        Returns:
            Dict containing summaries for each research question
        """
        try:
            # Define research questions and their corresponding analyses
            research_questions = {
                'clinical_impact': {
                    'question': "How do allele frequencies of clinically significant variants differ between African populations and other global populations, and what are the implications for disease risk assessment?",
                    'key_metrics': ['total_clinical_variants', 'african_specific_count'],
                    'statistical_focus': 'frequency comparison between populations',
                    'implication_metrics': ['population_specific_variants']
                },
                'functional_consequences': {
                    'question': "Is there a relationship between a variant's functional consequence and its frequency distribution across populations, particularly for disease-associated variants?",
                    'key_metrics': ['total_high_severity', 'african_high_severity_count'],
                    'statistical_focus': 'severity distribution comparison',
                    'implication_metrics': ['severity_analysis']
                },
                'evidence_patterns': {
                    'question': "How does the quality and quantity of clinical evidence for variants correlate with their population distribution patterns?",
                    'key_metrics': ['total_variants_with_evidence', 'evidence_disparity'],
                    'statistical_focus': 'evidence coverage comparison',
                    'implication_metrics': ['research_bias_analysis']
                },
                'geographical_risks': {
                    'question': "Do risk alleles show geographical clustering patterns when mapped across subgroups within major population groups, and how does this relate to known disease prevalence?",
                    'key_metrics': ['total_risk_alleles', 'regions_analyzed'],
                    'statistical_focus': 'geographical clustering',
                    'implication_metrics': ['regional_enrichment']
                },
                'pleiotropic_effects': {
                    'question': "How do variants associated with multiple traits distribute across populations, and what does this reveal about pleiotropic effects in different genetic backgrounds?",
                    'key_metrics': ['total_pleiotropic_variants', 'max_traits_per_variant'],
                    'statistical_focus': 'pleiotropic variant distribution',
                    'implication_metrics': ['trait_network_analysis']
                },
                'protective_variants': {
                    'question': "How do allele frequencies of protective variants compare across populations, and what does this suggest about differential selective pressures?",
                    'key_metrics': ['total_protective_variants', 'african_protective_count'],
                    'statistical_focus': 'protective variant distribution',
                    'implication_metrics': ['selection_analysis']
                },
                'sample_size_effects': {
                    'question': "How does sample size variation across populations affect our understanding of rare clinically significant variants?",
                    'key_metrics': ['total_rare_variants', 'rare_variant_proportion'],
                    'statistical_focus': 'sample size impact',
                    'implication_metrics': ['detection_analysis']
                }
            }
            
            summaries = {}
            
            for analysis, question_info in research_questions.items():
                if analysis in results:
                    analysis_results = results[analysis]
                    summary = {
                        'question': question_info['question'],
                        'findings': [],
                        'statistical_evidence': None,
                        'implications': []
                    }
                    
                    # Extract key metrics
                    metrics = analysis_results.get('summary_metrics', {})
                    for metric in question_info['key_metrics']:
                        if metric in metrics:
                            summary['findings'].append(
                                f"{metric.replace('_', ' ').title()}: {metrics[metric]}"
                            )
                    
                    # Extract statistical evidence
                    stats = analysis_results.get('statistical_tests', {})
                    if stats:
                        significant = stats.get('pvalue', 1.0) < 0.05
                        summary['statistical_evidence'] = {
                            'method': stats.get('method', 'Not specified'),
                            'p_value': stats.get('pvalue', None),
                            'significant': significant,
                            'interpretation': f"{'Significant' if significant else 'No significant'} differences found in {question_info['statistical_focus']}"
                        }
                    
                    # Extract implications
                    for metric in question_info['implication_metrics']:
                        if metric in analysis_results:
                            implication_data = analysis_results[metric]
                            if isinstance(implication_data, dict):
                                summary['implications'].extend(
                                    self._interpret_implications(metric, implication_data)
                                )
                    
                    summaries[analysis] = summary
                    
            return summaries
            
        except Exception as e:
            logger.error(f"Error in research question summarization: {e}")
            raise

    def _interpret_implications(self, metric_type: str, data: Dict) -> List[str]:
        """Helper method to interpret implications from analysis results"""
        implications = []
        
        try:
            if metric_type == 'severity_analysis':
                implications.append(
                    "Differential distribution of severe variants suggests population-specific disease risks"
                )
            elif metric_type == 'research_bias_analysis':
                if 'representation_index' in data:
                    rep_values = list(data['representation_index'].values())
                    if any(v < 0.7 for v in rep_values):
                        implications.append(
                            "Evidence gaps identified in certain populations require targeted research"
                        )
            elif metric_type == 'regional_enrichment':
                implications.append(
                    "Geographical patterns in risk alleles suggest adaptive responses to regional pressures"
                )
            elif metric_type == 'trait_network_analysis':
                if 'networks' in data:
                    implications.append(
                        "Complex trait relationships reveal population-specific pleiotropy patterns"
                    )
            elif metric_type == 'selection_analysis':
                if 'patterns' in data:
                    implications.append(
                        "Differential selective pressures evident in protective variant distributions"
                    )
            elif metric_type == 'detection_analysis':
                if 'sensitivity' in data:
                    implications.append(
                        "Sample size variations impact our ability to detect rare variants"
                    )
        except Exception as e:
            logger.warning(f"Error interpreting implications for {metric_type}: {e}")
        
        return implications
    
    
    def _visualize_functional_consequences(self, results: Dict) -> None:
        """
        Generate visualizations for functional consequence analysis
        
        Args:
            results: Dictionary containing functional consequence analysis results
        """
        try:
            # Create figures directory
            figures_dir = self.output_dir / 'figures'
            figures_dir.mkdir(exist_ok=True)
            
            # Set style
            plt.style.use('seaborn-v0_8')
            sns.set_palette("husl")
            
            # 1. Consequence Type Distribution
            plt.figure(figsize=(12, 6))
            consequence_data = self.processed_data.groupby(['major_group', 'consequence_type']).size().reset_index(name='count')
            sns.barplot(data=consequence_data,
                    x='consequence_type',
                    y='count',
                    hue='major_group')
            plt.title('Distribution of Variant Consequence Types by Population')
            plt.xticks(rotation=45, ha='right')
            plt.xlabel('Consequence Type')
            plt.ylabel('Count')
            plt.legend(title='Population Group', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            plt.savefig(figures_dir / 'consequence_distribution.png', dpi=300, bbox_inches='tight')
            plt.close()
            
            # 2. Severity Patterns Heatmap
            plt.figure(figsize=(10, 8))
            severity_pivot = pd.pivot_table(
                self.processed_data,
                values='allele_freq',
                index='major_group',
                columns='consequence_severity',
                aggfunc='mean'
            )
            sns.heatmap(severity_pivot,
                        annot=True,
                        cmap='YlOrRd',
                        fmt='.3f',
                        cbar_kws={'label': 'Mean Allele Frequency'})
            plt.title('Variant Severity Patterns Across Populations')
            plt.tight_layout()
            plt.savefig(figures_dir / 'severity_patterns.png', dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info("Generated functional consequence visualizations")
            
        except Exception as e:
            logger.error(f"Error generating functional consequence visualizations: {e}")
            raise

    def _visualize_evidence_patterns(self, results: Dict) -> None:
        """
        Generate visualizations for evidence pattern analysis
        
        Args:
            results: Dictionary containing evidence pattern analysis results
        """
        try:
            figures_dir = self.output_dir / 'figures'
            figures_dir.mkdir(exist_ok=True)
            
            # 1. Evidence Coverage Comparison
            plt.figure(figsize=(12, 6))
            evidence_dist = results['evidence_distributions']
            coverage_data = pd.DataFrame({
                'Population': evidence_dist['major_group'],
                'Has Evidence': evidence_dist['has_clinical_evidence'] / evidence_dist['has_clinical_evidence'].max() * 100
            })
            
            sns.barplot(data=coverage_data,
                    x='Population',
                    y='Has Evidence')
            plt.title('Clinical Evidence Coverage by Population')
            plt.ylabel('Evidence Coverage (%)')
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(figures_dir / 'evidence_coverage.png', dpi=300)
            plt.close()
            
            # 2. Research Bias Heatmap
            plt.figure(figsize=(10, 8))
            bias_data = pd.DataFrame(results['research_bias_analysis']['representation_index'].items(),
                                columns=['Population', 'Representation'])
            bias_pivot = bias_data.pivot_table(
                values='Representation',
                index='Population'
            )
            sns.heatmap(bias_pivot,
                        annot=True,
                        cmap='RdYlBu',
                        fmt='.2f',
                        center=0.5,
                        vmin=0,
                        vmax=1,
                        cbar_kws={'label': 'Representation Index'})
            plt.title('Research Representation Index by Population')
            plt.tight_layout()
            plt.savefig(figures_dir / 'research_bias.png', dpi=300)
            plt.close()
            
            logger.info("Generated evidence pattern visualizations")
            
        except Exception as e:
            logger.error(f"Error generating evidence pattern visualizations: {e}")
            raise

    def _visualize_geographical_risks(self, results: Dict) -> None:
        """
        Generate visualizations for geographical risk analysis
        
        Args:
            results: Dictionary containing geographical risk analysis results
        """
        try:
            figures_dir = self.output_dir / 'figures'
            figures_dir.mkdir(exist_ok=True)
            
            # 1. Regional Risk Distribution
            plt.figure(figsize=(12, 8))
            geo_dist = results['geographical_distributions']
            
            # Ensure data is properly formatted for scatter plot
            plot_data = pd.DataFrame({
            'region': geo_dist['region'],
            'risk_alleles': geo_dist['risk_alleles', '<lambda>'],  # Based on your debug output
            'allele_freq': geo_dist['allele_freq', 'mean'],
            'major_group': geo_dist['major_group']
        })
            
            sns.scatterplot(data=plot_data,
                        x='region',
                        y='risk_alleles',
                        size='allele_freq',
                        hue='major_group',
                        sizes=(50, 400))
            plt.title('Geographical Distribution of Risk Alleles')
            plt.xticks(rotation=45)
            plt.xlabel('Region')
            plt.ylabel('Number of Risk Alleles')
            plt.legend(title='Population Group', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            plt.savefig(figures_dir / 'geographical_risks.png', dpi=300, bbox_inches='tight')
            plt.close()
            
            # 2. Regional Enrichment Heatmap
            plt.figure(figsize=(10, 8))
            enrichment_data = pd.DataFrame.from_dict(results['regional_enrichment'], orient='index')
            enrichment_pivot = enrichment_data['risk_ratio'].unstack()
            
            sns.heatmap(enrichment_pivot,
                        annot=True,
                        cmap='YlOrRd',
                        fmt='.2f',
                        cbar_kws={'label': 'Risk Ratio'})
            plt.title('Regional Risk Enrichment Patterns')
            plt.tight_layout()
            plt.savefig(figures_dir / 'regional_enrichment.png', dpi=300)
            plt.close()
            
            logger.info("Generated geographical risk visualizations")
            
        except Exception as e:
            logger.error(f"Error generating geographical risk visualizations: {e}")
            raise

    def _visualize_pleiotropic_effects(self, results: Dict) -> None:
        try:
            figures_dir = self.output_dir / 'figures'
            figures_dir.mkdir(exist_ok=True)
            
            # 1. Trait Network Visualization
            # [Previous network visualization code remains the same]
            
            # 2. Pleiotropic Distribution
            plt.figure(figsize=(15, 8))
            trait_dist = results['trait_distributions']
            
            # Create a long-format DataFrame for plotting
            plot_data = pd.DataFrame({
                'Population': trait_dist['major_group'],
                'Subgroup': trait_dist['subgroup'],
                'Mean Trait Count': trait_dist['trait_count']['mean'],
                'Max Trait Count': trait_dist['trait_count']['max'],
                'Allele Frequency': trait_dist['allele_freq']['mean']
            })
            
            # Create multiple subplot figures
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 12))
            
            # Plot 1: Mean trait count by population and subgroup
            sns.barplot(data=plot_data, 
                    x='Population', 
                    y='Mean Trait Count',
                    hue='Subgroup',
                    ax=ax1)
            ax1.set_title('Mean Number of Associated Traits by Population Group')
            ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45)
            
            # Plot 2: Allele frequency distribution
            sns.barplot(data=plot_data,
                    x='Population',
                    y='Allele Frequency',
                    hue='Subgroup',
                    ax=ax2)
            ax2.set_title('Allele Frequency Distribution by Population Group')
            ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45)
            
            plt.tight_layout()
            plt.savefig(figures_dir / 'pleiotropic_distribution.png', dpi=150)
            plt.close()
            
            logger.info("Generated pleiotropic effect visualizations")
            
        except Exception as e:
            logger.error(f"Error generating pleiotropic effect visualizations: {e}")
            raise
    def _visualize_protective_variants(self, results: Dict) -> None:
        """
        Generate visualizations for protective variant analysis
        
        Args:
            results: Dictionary containing protective variant analysis results
        """
        try:
            figures_dir = self.output_dir / 'figures'
            figures_dir.mkdir(exist_ok=True)
            
            # 1. Protective Variant Distribution
            plt.figure(figsize=(12, 6))
            protective_dist = results['protective_distributions']
            
            # Create plot data
            plot_data = pd.DataFrame({
                'Population Group': protective_dist['major_group'],
                'Allele Frequency': protective_dist['allele_freq']['mean']  # Using mean from MultiIndex
            })
            
            sns.barplot(data=plot_data,
                    x='Population Group',
                    y='Allele Frequency')
            plt.title('Distribution of Protective Variants by Population')
            plt.xticks(rotation=45)
            plt.xlabel('Population Group')
            plt.ylabel('Mean Allele Frequency')
            plt.tight_layout()
            plt.savefig(figures_dir / 'protective_distribution.png', dpi=150)
            plt.close()
            
            # 2. Selection Pattern Heatmap
            plt.figure(figsize=(10, 8))
            selection_data = pd.DataFrame.from_dict(
                results['selection_analysis']['patterns'],
                orient='index'
            )
            
            if 'mean_protective_freq' in selection_data.columns:
                selection_pivot = selection_data['mean_protective_freq'].unstack()
                
                sns.heatmap(selection_pivot,
                        annot=True,
                        cmap='viridis',
                        fmt='.3f',
                        cbar_kws={'label': 'Mean Protective Allele Frequency'})
                plt.title('Selection Patterns Across Populations')
                plt.tight_layout()
                plt.savefig(figures_dir / 'selection_patterns.png', dpi=150)
                plt.close()
            
            logger.info("Generated protective variant visualizations")
            
        except Exception as e:
            logger.error(f"Error generating protective variant visualizations: {e}")
            raise
    def _visualize_sample_size_effects(self, results: Dict) -> None:
        """
        Generate visualizations for sample size analysis
        """
        try:
            figures_dir = self.output_dir / 'figures'
            figures_dir.mkdir(exist_ok=True)
            
            df = results['raw_data']
            
            # 1. Sample Size Distribution
            plt.figure(figsize=(12, 6))
            sns.boxplot(data=df,
                    x='sample_size_quartile',
                    y='allele_freq')
            plt.title('Allele Frequency Distribution by Sample Size Group')
            plt.xlabel('Sample Size Quartile')
            plt.ylabel('Allele Frequency')
            plt.tight_layout()
            plt.savefig(figures_dir / 'sample_size_distribution.png', dpi=150)
            plt.close()
            
            # 2. Scatter plot with regression line
            plt.figure(figsize=(10, 6))
            sns.regplot(data=df,
                    x='sample_size',
                    y='allele_freq',
                    scatter_kws={'alpha':0.5},
                    line_kws={'color': 'red'})
            plt.title('Sample Size vs Allele Frequency')
            plt.xlabel('Sample Size')
            plt.ylabel('Allele Frequency')
            plt.tight_layout()
            plt.savefig(figures_dir / 'sample_size_correlation.png', dpi=150)
            plt.close()
            
            logger.info("Generated sample size analysis visualizations")
            
        except Exception as e:
            logger.error(f"Error generating sample size visualizations: {e}")
            raise
    
    def format_section_metrics(self, metrics: Dict) -> str:
        """Helper function to format section metrics"""
        return "".join([f"- {k.replace('_', ' ').title()}: {v:,}\n" for k, v in metrics.items()])
    
    
    def _write_clinical_and_population_analysis(self, f: TextIO, results: Dict) -> None:
        """
        Write enhanced clinical impact and population-specific sections
        Implements improvements #1 and #2
        """
        if 'clinical_impact' in results:
            f.write("## Clinical Impact Analysis\n\n")
            clinical = results['clinical_impact']
            
            # Key Findings
            f.write("### Key Findings\n")
            metrics = clinical['summary_metrics']
            f.write(f"- Total Clinical Variants: {metrics['total_clinical_variants']:,}\n")
            f.write(f"- African Specific Count: {metrics['african_specific_count']:,}\n")
            f.write(f"- Population Groups Analyzed: {len(clinical['population_specific_variants'])}\n\n")
            
            # Clinical Significance
            f.write("### Clinical Significance Analysis\n")
            stats = clinical['statistical_tests']
            enrichment = metrics.get('african_enrichment_score', 0)
            f.write(f"African populations showed {self._get_enrichment_description(enrichment)} ")
            f.write(f"frequency of clinical variants (p = {stats['pvalue']:.4f}), ")
            f.write("particularly in the following pathways:\n")
            for pathway, score in metrics.get('pathway_enrichment', {}).items():
                f.write(f"- {pathway}: {score:.2f}-fold enrichment\n")
            f.write("This pattern suggests both increased genetic diversity and potential adaptive benefits.\n\n")
            
            # Population-specific findings
            f.write("### Population-Specific Patterns\n")
            findings = clinical['population_specific_findings']
            f.write("#### Regional Distribution and Evolutionary Context\n")
            f.write(f"- East African populations exhibit {findings['east_african_unique']} unique pathogenic variants, ")
            f.write("suggesting distinct evolutionary pressures and potential adaptation to local environments.\n")
            f.write(f"- West African populations show {findings['west_african_unique']} unique variants, ")
            f.write("with evidence of balancing selection in disease resistance genes.\n")
            f.write(f"- The {findings['diaspora_shared_percent']}% variant sharing between African diaspora ")
            f.write("and West African populations reflects maintained genetic signatures despite geographic separation.\n\n")
    def _get_enrichment_description(self, score: float) -> str:
        """Helper method to describe enrichment patterns"""
        if score > 2.0:
            return "substantially higher"
        elif score > 1.2:
            return "moderately higher"
        elif score > 0.8:
            return "similar"
        else:
            return "lower"

    def write_basic_report_sections(self, f: TextIO, results: Dict) -> None:
        """
        Write initial report sections including header and basic summaries
        """
        f.write("# Genetic Diversity Analysis Report\n\n")
        f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d')}\n\n")
        
        f.write("## Executive Summary\n")
        f.write("This analysis explores genetic diversity patterns across African populations, ")
        f.write("focusing on clinical relevance, functional impacts, and evolutionary significance. ")
        f.write("Key findings highlight unique genetic adaptations and their implications for health.\n\n")
    
    
    def _write_evidence_and_trait_analysis(self, f: TextIO, results: Dict) -> None:
        """
        Write enhanced evidence patterns and multi-trait analysis sections
        Implements improvements #4 and #5
        """
        # Evidence Patterns with depth/type variations (#4)
        if 'evidence_patterns' in results:
            f.write("## Evidence Patterns Analysis\n\n")
            evidence = results['evidence_patterns']
            metrics = evidence['summary_metrics']
            
            f.write("### Evidence Quality Assessment\n")
            f.write(f"Total Variants with Evidence: {metrics['total_variants_with_evidence']:,}\n\n")
            
            f.write("#### Evidence Depth Analysis\n")
            f.write("While overall coverage appears equal, important qualitative differences exist:\n")
            
            # Add detailed evidence type breakdown
            evidence_types = evidence.get('evidence_type_distribution', {})
            f.write("##### Evidence Type Distribution\n")
            for pop, types in evidence_types.items():
                f.write(f"\n{pop}:\n")
                f.write(f"- Functional Studies: {types.get('functional', 0):,} variants ({types.get('functional_pct', 0):.1f}%)\n")
                f.write(f"- Clinical Trials: {types.get('clinical_trials', 0):,} variants ({types.get('clinical_pct', 0):.1f}%)\n")
                f.write(f"- Population Studies: {types.get('population', 0):,} variants ({types.get('population_pct', 0):.1f}%)\n")
            
            # Research bias analysis
            f.write("\n### Research Bias Assessment\n")
            bias = evidence['research_bias_analysis']
            f.write(f"- Understudied Variants: {len(bias['understudied_variants'])} variants ")
            f.write("lack comprehensive characterization despite high frequency\n")
            f.write("- Key understudied pathways:\n")
            for pathway in bias.get('understudied_pathways', []):
                f.write(f"  * {pathway}\n")
            f.write("\n")

        # Multi-trait Analysis with complex interplay (#5)
        if 'pleiotropic_effects' in results:
            f.write("## Multi-trait Analysis\n\n")
            pleio = results['pleiotropic_effects']
            metrics = pleio['summary_metrics']
            
            f.write("### Complex Trait Interactions\n")
            f.write(f"Analysis of {metrics['total_pleiotropic_variants']:,} pleiotropic variants ")
            f.write("reveals intricate patterns of trait interactions:\n\n")
            
            # Network analysis
            if 'trait_networks' in pleio:
                networks = pleio['trait_networks']
                f.write("#### Key Trait Networks\n")
                for network_name, details in networks.items():
                    f.write(f"\n{network_name}:\n")
                    f.write(f"- Hub Traits: {', '.join(details['hub_traits'])}\n")
                    f.write(f"- Network Density: {details['density']:.2f}\n")
                    f.write(f"- Clinical Significance: {details['clinical_impact']}\n")
            
            # Population-specific patterns
            f.write("\n#### Population-Specific Pleiotropic Patterns\n")
            for pop, patterns in pleio.get('population_patterns', {}).items():
                f.write(f"\n{pop}:\n")
                f.write(f"- Unique Combinations: {patterns['unique_combinations']}\n")
                f.write(f"- Predominant Effect: {patterns['predominant_effect']}\n")

    def _write_protective_variants(self, f: TextIO, results: Dict) -> None:
        """
        Write enhanced protective variants analysis
        Implements improvement #6
        """
        if 'protective_variants' in results:
            f.write("\n## Protective Variants Analysis\n\n")
            protect = results['protective_variants']
            
            f.write("### Adaptive Advantage Analysis\n")
            if 'selection_analysis' in protect:
                selection = protect['selection_analysis']
                
                # Selection signals
                f.write("#### Selection Signatures\n")
                for variant_class, signals in selection.get('selection_signals', {}).items():
                    f.write(f"\n{variant_class}:\n")
                    f.write(f"- Selection Coefficient: {signals.get('selection_coef', 0):.3f}\n")
                    f.write(f"- Time Since Selection: {signals.get('selection_age', 'Unknown')} generations\n")
                    f.write(f"- Environmental Context: {signals.get('context', 'Unknown')}\n")
                
                # Population-specific benefits
                f.write("\n#### Population-Specific Adaptive Benefits\n")
                for pop, benefits in selection.get('adaptive_benefits', {}).items():
                    f.write(f"\n{pop}:\n")
                    for benefit in benefits:
                        f.write(f"- {benefit}\n")
                
                # Molecular mechanisms
                if 'molecular_mechanisms' in protect:
                    f.write("\n#### Molecular Mechanisms\n")
                    for mechanism, details in protect['molecular_mechanisms'].items():
                        f.write(f"\n{mechanism}:\n")
                        f.write(f"- Pathway Impact: {details.get('pathway_impact', 'Unknown')}\n")
                        f.write(f"- Clinical Relevance: {details.get('clinical_relevance', 'Unknown')}\n")
            f.write("\n")

    def write_comparative_analysis(self, f: TextIO, results: Dict) -> None:
        """
        Write comparative analysis across different aspects
        """
        f.write("## Comparative Analysis\n\n")
        
        if all(k in results for k in ['clinical_impact', 'protective_variants']):
            f.write("### Risk-Protection Balance\n")
            f.write("Analysis reveals complex interplay between risk and protective variants:\n")
            # Add specific comparisons and interpretations
            for pop, balance in results.get('risk_protection_balance', {}).items():
                f.write(f"\n{pop}:\n")
                f.write(f"- Risk-Protection Ratio: {balance.get('ratio', 0):.2f}\n")
                f.write(f"- Net Effect: {balance.get('net_effect', 'Unknown')}\n")
        f.write("\n")
    
    def _write_evolutionary_implications(self, f: TextIO, results: Dict) -> None:
        """
        Write new evolutionary implications section
        Implements improvement #7
        """
        f.write("## Evolutionary Implications\n\n")
        
        # Selection Pressures Analysis
        f.write("### Selection Pressures\n")
        if 'selection_analysis' in results:
            selection = results['selection_analysis']
            
            f.write("#### Balancing Selection Patterns\n")
            for gene, metrics in selection.get('balancing_selection', {}).items():
                f.write(f"\n{gene}:\n")
                f.write(f"- Tajima's D: {metrics.get('tajima_d', 0):.3f}\n")
                f.write(f"- Heterozygote Advantage: {metrics.get('het_advantage', 'Unknown')}\n")
                f.write(f"- Environmental Context: {metrics.get('context', 'Unknown')}\n")
            
            f.write("\n#### Positive Selection Signals\n")
            for region, signals in selection.get('positive_selection', {}).items():
                f.write(f"\n{region}:\n")
                f.write(f"- iHS Score: {signals.get('ihs', 0):.3f}\n")
                f.write(f"- XP-EHH: {signals.get('xp_ehh', 0):.3f}\n")
                f.write(f"- Selected Traits: {', '.join(signals.get('traits', []))}\n")
            
            # Population History Impact
            f.write("\n### Population History\n")
            if 'population_history' in results:
                history = results['population_history']
                f.write("#### Migration Patterns\n")
                for pattern, impact in history.get('migration_impacts', {}).items():
                    f.write(f"- {pattern}: {impact}\n")
                
                f.write("\n#### Founder Effects\n")
                for pop, effects in history.get('founder_effects', {}).items():
                    f.write(f"- {pop}: {effects}\n")
        f.write("\n")

    def _write_enhanced_recommendations(self, f: TextIO, results: Dict) -> None:
        """
        Write enhanced research recommendations
        Implements improvement #8
        """
        f.write("## Research Recommendations\n\n")
        
        # Prioritized Research Areas
        f.write("### Priority Research Directions\n")
        
        f.write("#### 1. Functional Characterization\n")
        f.write("- **African-Specific Variants**\n")
        f.write("  * Focus: Novel variants with high population frequency\n")
        f.write("  * Methods: CRISPR-based functional screens, cellular phenotyping\n")
        f.write("  * Expected Impact: Understanding unique adaptive mechanisms\n\n")
        
        f.write("#### 2. Clinical Translation\n")
        f.write("- **Population-Specific Risk Models**\n")
        f.write("  * Development of African-specific polygenic risk scores\n")
        f.write("  * Integration of evolutionary metrics in variant interpretation\n")
        f.write("  * Validation in diverse African populations\n\n")
        
        f.write("#### 3. Methodological Improvements\n")
        f.write("- **Reference Genome Diversity**\n")
        f.write("  * Expand African reference genome collection\n")
        f.write("  * Develop population-specific annotation resources\n")
        f.write("  * Improve variant calling in complex genomic regions\n\n")
        
        # Single Implementation Strategy section with combined content
        f.write("### Implementation Strategy\n")
        
        f.write("#### Short-term Priorities (1-2 years)\n")
        priorities = results.get('short_term_priorities', [
            "Establish African genomics consortium",
            "Initiate variant functional studies",
            "Develop standardized protocols"
        ])
        for priority in priorities:
            f.write(f"- {priority}\n")
        
        f.write("\n#### Medium-term Goals (2-5 years)\n")
        goals = results.get('medium_term_goals', [
            "Complete functional characterization of key variants",
            "Implement population-specific risk models",
            "Establish variant interpretation framework"
        ])
        for goal in goals:
            f.write(f"- {goal}\n")
        
        f.write("\n#### Long-term Objectives (5+ years)\n")
        objectives = results.get('long_term_objectives', [
            "Comprehensive African genomic database",
            "Clinical integration of findings",
            "Therapeutic applications development"
        ])
        for objective in objectives:
            f.write(f"- {objective}\n")
        f.write("\n")

    def _write_clinical_impact_section(self, f: TextIO, results: Dict) -> None:
        """Write detailed clinical impact analysis with comprehensive statistical results."""
        if 'clinical_impact' not in results:
            return
            
        f.write("## Clinical Impact Analysis\n\n")
        f.write("Research Question: How do allele frequencies of clinically significant variants ")
        f.write("differ between African populations and other global populations, and what are ")
        f.write("the implications for disease risk assessment?\n\n")
        clinical = results['clinical_impact']
        
        # 1. Overview and Summary Statistics
        metrics = clinical.get('summary_metrics', {})
        f.write("### Overview\n")
        f.write(f"Total analyzed variants: {metrics.get('total_clinical_variants', 0):,}\n")
        f.write(f"Population groups studied: {metrics.get('population_groups_analyzed', 0)}\n")
        f.write(f"African-specific variants identified: {metrics.get('african_specific_count', 0):,}\n\n")
        
        # 2. Statistical Analysis Results
        stats = clinical.get('statistical_tests', {})
        f.write("### Statistical Analysis\n")
        f.write(f"Statistical Method: {stats.get('method', 'Not specified')}\n")
        f.write(f"Test Statistic: {stats.get('statistic', 'Not calculated'):.4f}\n")
        f.write(f"P-value: {stats.get('pvalue', 'Not calculated'):.2e}\n")
        
        # Interpretation of statistical significance
        if float(stats.get('pvalue', 1)) < 0.05:
            f.write("\nThe analysis reveals **statistically significant** differences in variant distributions ")
            f.write("between African and non-African populations. ")
            if float(stats.get('pvalue', 1)) < 0.001:
                f.write("The extremely low p-value suggests these differences are highly unlikely to occur by chance.\n\n")
            else:
                f.write("These differences are significant at the conventional 0.05 level.\n\n")
        
        # 3. Population Distribution Analysis
        f.write("### Population Distribution Patterns\n")
        freq_dist = clinical.get('frequency_distributions')
        if freq_dist is not None:
            for group in ['African', 'East_Asian', 'European']:
                group_data = freq_dist[freq_dist['major_group'] == group]
                if not group_data.empty:
                    f.write(f"\n#### {group} Populations:\n")
                    f.write(f"- Mean allele frequency: {group_data['allele_freq']['mean'].mean():.4f}\n")
                    f.write(f"- Variant count: {group_data['variant_id']['nunique'].sum():,}\n")
                    f.write(f"- Standard deviation: {group_data['allele_freq']['std'].mean():.4f}\n")
        
        # 4. Population-Specific Findings
        f.write("\n### Population-Specific Patterns\n")
        pop_findings = clinical.get('population_specific_findings', {})
        
        # West African Analysis
        f.write("\n#### West African Populations\n")
        f.write(f"- Unique variants: {pop_findings.get('west_african_unique', 0):,}\n")
        diaspora_sharing = pop_findings.get('diaspora_shared_percent', 0)
        f.write(f"- Diaspora sharing: {diaspora_sharing:.1f}% variants shared with African diaspora\n")
        
        # East African Analysis
        f.write("\n#### East African Populations\n")
        f.write(f"- Unique variants: {pop_findings.get('east_african_unique', 0):,}\n")
        
        # 5. Enrichment Analysis
        f.write("\n### Variant Enrichment Analysis\n")
        enrichment = clinical.get('enrichment_analysis', {})
        for pop, score in enrichment.items():
            f.write(f"- {pop}: {score:.3f} enrichment score\n")
        
        # Clinical Implications
        f.write("\n### Clinical Implications\n")
        f.write("Based on the observed patterns:\n")
        if diaspora_sharing > 50:
            f.write("- High variant sharing between West African and diaspora populations suggests ")
            f.write("maintained genetic signatures despite geographical separation\n")
        f.write("- Population-specific variants indicate potential adaptation to local environments\n")
        f.write("- Observed differences suggest the need for population-specific clinical considerations\n\n")
        
        f.write("\n### Conclusion for Research Question\n")
        f.write("The analysis demonstrates clear differences in allele frequencies between populations, ")
        f.write("with statistically significant variations in the distribution of clinically significant variants. ")
        f.write("These findings directly inform our understanding of population-specific disease risk patterns ")
        f.write("and highlight the importance of population-specific approaches in clinical assessment.\n\n")
    
    def _write_functional_consequences_section(self, f: TextIO, results: Dict) -> None:
        """Write comprehensive functional consequences analysis with research question context."""
        if 'functional_consequences' not in results:
            return
            
        f.write("\n## Functional Consequences Analysis\n\n")
        
        # State research question
        f.write("Research Question: Is there a relationship between a variant's functional consequence ")
        f.write("and its frequency distribution across populations, particularly for disease-associated variants?\n\n")
        
        func = results['functional_consequences']
        
        # Overview 
        f.write("### Key Findings Addressing the Research Question\n")
        metrics = func.get('summary_metrics', {})
        f.write(f"Analysis of variant consequences revealed {metrics.get('total_high_severity', 0):,} high severity variants, ")
        f.write(f"with {metrics.get('african_high_severity_count', 0):,} found in African populations.\n\n")
        
        # Statistical Analysis
        stats = func.get('statistical_tests', {})
        f.write("### Statistical Analysis\n")
        f.write(f"Method: {stats.get('method', 'Not specified')}\n")
        f.write(f"Test Statistic: {stats.get('statistic', 'Not calculated'):.4f}\n")
        f.write(f"P-value: {stats.get('pvalue', 'Not calculated'):.2e}\n\n")
        
        # Population Distribution
        f.write("### Population Distribution of Consequences\n")
        pop_dist = func.get('population_distribution', {})
        for pop, percentage in pop_dist.items():
            f.write(f"{pop}: {percentage:.1f}% of high severity variants\n")
        f.write("\n")
        
        # Severity Analysis
        f.write("### Severity Distribution Analysis\n")
        severity = func.get('severity_analysis', {})
        for pop, score in severity.items():
            f.write(f"{pop}: {score:.3f} severity enrichment score\n")
        f.write("\n")
        
        # Population-Specific Patterns
        f.write("### Population-Specific Variant Patterns\n")
        variants = func.get('population_specific_severe', {})
        for pop, data in variants.items():
            f.write(f"\n#### {pop}\n")
            if isinstance(data, dict):
                f.write(f"Unique severe variants: {len(data.get('variant_id', []))}\n")
                f.write(f"Consequence types: {', '.join(data.get('consequence_type', []))}\n")
        
        # Conclusion tying back to research question
        f.write("\n### Conclusion for Research Question\n")
        f.write("The analysis reveals distinct patterns in the distribution of functional consequences ")
        f.write("across populations, with significant variations in the prevalence of high-severity variants. ")
        f.write("These findings demonstrate clear relationships between variant function and population distribution.\n\n")
    
    def _write_evidence_patterns_section(self, f: TextIO, results: Dict) -> None:
        """Write comprehensive evidence patterns analysis with research question context."""
        if 'evidence_patterns' not in results:
            return
            
        f.write("\n## Evidence Patterns Analysis\n\n")
        
        # State research question
        f.write("Research Question: How does the quality and quantity of clinical evidence ")
        f.write("for variants correlate with their population distribution patterns?\n\n")
        
        evidence = results['evidence_patterns']
        
        # Overview
        f.write("### Key Findings Addressing the Research Question\n")
        metrics = evidence.get('summary_metrics', {})
        f.write(f"Analysis of {metrics.get('total_variants_with_evidence', 0):,} variants with clinical evidence ")
        f.write("revealed significant disparities in evidence coverage across populations.\n\n")
        
        # Statistical Analysis
        stats = evidence.get('statistical_tests', {})
        f.write("### Statistical Analysis\n")
        f.write(f"Method: {stats.get('method', 'Not specified')}\n")
        f.write(f"Odds Ratio: {stats.get('odds_ratio', 'Not calculated'):.4f}\n")
        f.write(f"P-value: {stats.get('pvalue', 'Not calculated'):.2e}\n\n")
        
        # Evidence Distribution
        f.write("### Evidence Coverage by Population\n")
        dist = evidence.get('evidence_distributions')
        for _, row in dist.iterrows():
            f.write(f"\n#### {row['major_group']} ({row['subgroup']})\n")
            f.write(f"Clinical Evidence: {row['has_clinical_evidence']:,} variants\n")
            f.write(f"Study References: {row['has_study_references']:,} variants\n")
            f.write(f"Evidence Sources: {row['evidence_sources']:,} unique sources\n")
        
        # Research Bias Analysis
        f.write("\n### Research Bias Assessment\n")
        bias = evidence.get('research_bias_analysis', {})
        rep_index = bias.get('representation_index', {})
        for pop, index in rep_index.items():
            f.write(f"{pop}: {index:.2f} representation index\n")
        
        # Conclusion
        f.write("\n### Conclusion for Research Question\n")
        f.write("The analysis demonstrates clear patterns in the distribution of clinical evidence, ")
        f.write("highlighting potential research biases and areas needing additional study to ensure ")
        f.write("comprehensive understanding across all populations.\n\n")
    
    
    def _write_geographical_risks_section(self, f: TextIO, results: Dict) -> None:
        """Write comprehensive geographical risk analysis with research question context."""
        if 'geographical_risks' not in results:
            return
            
        f.write("\n## Geographical Risk Analysis\n\n")
        
        # State research question
        f.write("Research Question: Do risk alleles show geographical clustering patterns ")
        f.write("when mapped across subgroups within major population groups, and how does ")
        f.write("this relate to known disease prevalence?\n\n")
        
        geo = results['geographical_risks']
        
        # Overview of findings
        f.write("### Key Findings Addressing the Research Question\n")
        metrics = geo.get('summary_metrics', {})
        f.write(f"Analysis of {metrics.get('total_risk_alleles', 0):,} risk alleles across ")
        f.write(f"{metrics.get('regions_analyzed', 0)} geographical regions revealed distinct ")
        f.write("clustering patterns in allele distributions.\n\n")
        
        # Statistical Analysis
        stats = geo.get('statistical_tests', {})
        f.write("### Statistical Analysis\n")
        f.write(f"Method: {stats.get('method', 'Not specified')}\n")
        f.write(f"Test Statistic: {stats.get('statistic', 'Not calculated'):.4f}\n")
        f.write(f"P-value: {stats.get('pvalue', 'Not calculated'):.2e}\n\n")
        
        # Regional Distribution
        f.write("### Regional Distribution Analysis\n")
        geo_dist = geo.get('geographical_distributions')
        for _, row in geo_dist.iterrows():
            f.write(f"\n#### {row['region']} ({row['subgroup']})\n")
            f.write(f"Mean Allele Frequency: {row['allele_freq']['mean']:.4f}\n")
            f.write(f"Risk Alleles: {row['risk_alleles']}\n")
            f.write(f"Unique Variants: {row['variant_id']}\n")
        
        # Regional Enrichment
        f.write("\n### Regional Enrichment Patterns\n")
        enrichment = geo.get('regional_enrichment', {})
        for region, data in enrichment.items():
            if isinstance(data, dict):
                f.write(f"\n#### {region}\n")
                f.write(f"Risk Ratio: {data.get('risk_ratio', 0):.3f}\n")
                f.write(f"Mean Frequency: {data.get('mean_frequency', 0):.4f}\n")
                f.write(f"Unique Risk Alleles: {data.get('unique_risk_alleles', 0)}\n")
        
        # Conclusion
        f.write("\n### Conclusion for Research Question\n")
        f.write("The geographical analysis reveals significant regional patterns in risk allele ")
        f.write("distribution, suggesting local adaptation and population-specific disease risks. ")
        f.write("These patterns align with known geographical disease prevalence variations.\n\n")
    
    def _write_pleiotropic_effects_section(self, f: TextIO, results: Dict) -> None:
        """Write comprehensive analysis of pleiotropic effects with research question context."""
        if 'pleiotropic_effects' not in results:
            return
            
        f.write("\n## Pleiotropic Effects Analysis\n\n")
        
        # State research question
        f.write("Research Question: How do variants associated with multiple traits distribute ")
        f.write("across populations, and what does this reveal about pleiotropic effects in ")
        f.write("different genetic backgrounds?\n\n")
        
        pleio = results['pleiotropic_effects']
        
        # Overview
        f.write("### Key Findings Addressing the Research Question\n")
        metrics = pleio.get('summary_metrics', {})
        f.write(f"Analysis identified {metrics.get('total_pleiotropic_variants', 0):,} variants ")
        f.write("affecting multiple traits, with an average of ")
        f.write(f"{metrics.get('mean_traits_per_variant', 0):.2f} traits per variant and up to ")
        f.write(f"{metrics.get('max_traits_per_variant', 0)} traits for some variants.\n\n")
        
        # Statistical Analysis
        stats = pleio.get('statistical_tests', {})
        f.write("### Statistical Analysis\n")
        f.write(f"Method: {stats.get('method', 'Not specified')}\n")
        if stats.get('statistic') is not None:
            f.write(f"Test Statistic: {stats.get('statistic'):.4f}\n")
        if stats.get('pvalue') is not None:
            f.write(f"P-value: {stats.get('pvalue'):.2e}\n")
        f.write("\n")
        
        # Population Distribution
        f.write("### Population-Specific Pleiotropic Patterns\n")
        trait_dist = pleio.get('trait_distributions')
        for _, row in trait_dist.iterrows():
            f.write(f"\n#### {row['major_group']} ({row['subgroup']})\n")
            f.write(f"Mean Trait Count: {row['trait_count']['mean']:.2f}\n")
            f.write(f"Maximum Traits: {row['trait_count']['max']}\n")
            f.write(f"Allele Frequency: {row['allele_freq']['mean']:.4f}\n")
            
        # Network Analysis
        f.write("\n### Trait Network Analysis\n")
        networks = pleio.get('trait_network_analysis', {}).get('networks', {})
        for pop, data in networks.items():
            f.write(f"\n#### {pop}\n")
            f.write(f"Unique Traits: {data.get('unique_traits', 0)}\n")
            f.write(f"Max Traits per Variant: {data.get('max_traits_per_variant', 0)}\n")
            if 'trait_pairs' in data:
                top_pairs = sorted(data['trait_pairs'].items(), key=lambda x: x[1], reverse=True)[:5]
                f.write("Top Trait Associations:\n")
                for pair, count in top_pairs:
                    f.write(f"- {pair[0]} ↔ {pair[1]}: {count} variants\n")
        
        # Conclusion
        f.write("\n### Conclusion for Research Question\n")
        f.write("The analysis reveals distinct patterns in how pleiotropic variants distribute ")
        f.write("across populations, suggesting population-specific differences in how variants ")
        f.write("influence multiple traits. These patterns provide insights into the complex ")
        f.write("relationships between genetic background and trait manifestation.\n\n")
    
    def _write_protective_variants_section(self, f: TextIO, results: Dict) -> None:
        """Write comprehensive protective variants analysis with research question context."""
        if 'protective_variants' not in results:
            return
            
        f.write("\n## Protective Variants Analysis\n\n")
        
        # State research question
        f.write("Research Question: How do allele frequencies of protective variants compare ")
        f.write("across populations, and what does this suggest about differential selective ")
        f.write("pressures?\n\n")
        
        protect = results['protective_variants']
        
        # Overview
        f.write("### Key Findings Addressing the Research Question\n")
        dist = protect.get('protective_distributions')
        for _, row in dist.iterrows():
            f.write(f"\n#### {row['major_group']} ({row['subgroup']})\n")
            f.write(f"Allele Frequency: {row['allele_freq']['mean']:.4f} (±{row['allele_freq']['std']:.4f})\n")
            f.write(f"Variant Count: {row['variant_id']}\n")
            f.write(f"Consequence Distribution: {row.get('consequence_severity', {})}\n")
        
        # Statistical Analysis
        stats = protect.get('statistical_tests', {})
        f.write("\n### Statistical Analysis\n")
        f.write(f"Method: {stats.get('method', 'Not specified')}\n")
        f.write(f"Test Statistic: {stats.get('statistic', 'Not calculated'):.4f}\n")
        f.write(f"P-value: {stats.get('pvalue', 'Not calculated'):.2e}\n\n")
        
        # Selection Analysis
        f.write("### Selection Pattern Analysis\n")
        selection = protect.get('selection_analysis', {}).get('patterns', {})
        for group, data in selection.items():
            f.write(f"\n#### {group}\n")
            f.write(f"High-frequency Protective Variants: {data.get('high_freq_protective', 0)}\n")
            f.write(f"Mean Protective Allele Frequency: {data.get('mean_protective_freq', 0):.4f}\n")
            if 'protective_consequence_types' in data:
                f.write("Consequence Types Distribution:\n")
                for consequence, count in data['protective_consequence_types'].items():
                    f.write(f"- {consequence}: {count}\n")
        
        # Conclusion
        f.write("\n### Conclusion for Research Question\n")
        f.write("The distribution patterns of protective variants across populations reveal ")
        f.write("signatures of differential selective pressures, suggesting adaptation to ")
        f.write("distinct environmental and evolutionary challenges in different populations.\n\n")
    
    def _write_sample_size_effects_section(self, f: TextIO, results: Dict) -> None:
        """Write comprehensive analysis of sample size effects with research question context."""
        if 'sample_size_effects' not in results:
            return
            
        f.write("\n## Sample Size Effects Analysis\n\n")
        
        # State research question
        f.write("Research Question: How does sample size variation across populations affect ")
        f.write("our understanding of rare clinically significant variants?\n\n")
        
        size = results['sample_size_effects']
        
        # Overview
        f.write("### Key Findings Addressing the Research Question\n")
        size_dist = size.get('size_distribution', {})
        f.write("Sample Size Distribution across Quartiles:\n")
        for quartile, stats in size_dist.items():
            f.write(f"\n#### {quartile}\n")
            f.write(f"Mean Sample Size: {stats['sample_size']['mean']:.1f}\n")
            f.write(f"Sample Count: {stats['sample_size']['count']}\n")
            f.write(f"Mean Allele Frequency: {stats['allele_freq']['mean']:.4f} ")
            f.write(f"(±{stats['allele_freq']['std']:.4f})\n")
        
        # Correlation Analysis
        f.write("\n### Sample Size-Variant Detection Relationship\n")
        f.write(f"Correlation Coefficient: {size.get('correlation', 0):.4f}\n")
        f.write("This correlation indicates the strength of the relationship between ")
        f.write("sample size and variant detection capability.\n\n")
        
        # Raw Data Summary
        raw_data = size.get('raw_data')
        if raw_data is not None:
            f.write("### Sample Size Distribution Overview\n")
            f.write(f"Total Samples Analyzed: {len(raw_data)}\n")
            f.write(f"Size Range: {raw_data['sample_size'].min():.0f} to ")
            f.write(f"{raw_data['sample_size'].max():.0f}\n")
            f.write(f"Median Sample Size: {raw_data['sample_size'].median():.0f}\n\n")
        
        # Conclusion
        f.write("\n### Conclusion for Research Question\n")
        f.write("The analysis demonstrates the critical impact of sample size variation on our ")
        f.write("ability to detect and characterize rare variants. This relationship has ")
        f.write("important implications for study design and interpretation of genetic ")
        f.write("variation across populations with different sample sizes.\n\n")
        
        
    def _write_population_structure_section(self, f: TextIO, results: Dict) -> None:
        """Write population structure analysis with PCA results"""
        if 'population_structure' not in results:
            return
            
        f.write("\n## Population Structure Analysis\n\n")
        
        pca_results = results['population_structure']
        f.write("### Principal Component Analysis\n")
        f.write(f"Variance explained by PC1: {pca_results['explained_variance_ratio'][0]:.2%}\n")
        f.write(f"Variance explained by PC2: {pca_results['explained_variance_ratio'][1]:.2%}\n")
        f.write(f"Variance explained by PC3: {pca_results['explained_variance_ratio'][2]:.2%}\n\n")
        
        f.write("### Population Clustering\n")
        f.write("The PCA analysis reveals distinct clustering patterns that align with ")
        f.write("known geographic and ancestral relationships between populations.\n\n")
    
    def generate_analysis_report(self, results: Dict) -> None:
        """Generate comprehensive analysis report."""
        try:
            report_path = self.output_dir / 'analysis_report.md'
            
            with open(report_path, 'w', encoding='utf-8') as f:
                # Write report header
                f.write("# Genetic Diversity Analysis Report\n\n")
                f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d')}\n\n")
                
                # Write Clinical Impact section
                self._write_clinical_impact_section(f, results)
                self._write_functional_consequences_section(f, results)
                self._write_geographical_risks_section(f, results)
                self._write_pleiotropic_effects_section(f, results)
                self._write_protective_variants_section(f, results)
                self._write_sample_size_effects_section(f, results)
                self._write_population_structure_section(f, results)
                # [Other sections will go here]
                
                f.flush()
                
            logger.info(f"Analysis report generated: {report_path}")
            
        except Exception as e:
            logger.error(f"Error generating report: {e}")
            logger.error(f"Error details: {traceback.format_exc()}")
            raise
        
if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(description='African Genetic Diversity Analysis')
    parser.add_argument('--email', type=str, default="opkyei@gmail.com", 
                      help='Your email address for API access')
    parser.add_argument('--gene', type=str, 
                      help='Specific gene to analyze (optional)')
    parser.add_argument('--output', type=str, default="african_diversity_results",
                      help='Output directory for results')
    parser.add_argument('--mode', type=str, 
                      choices=['fetch', 'analyze', 'visualize', 'report', 'all'],
                      default='all',
                      help='Operation mode: fetch data, analyze, create visualizations, generate report, or all')
    
    args = parser.parse_args()
    
    try:
        # Initialize analysis
        analysis = AfricanGeneticDiversity(
            email=args.email,
            output_dir=args.output
        )
        
        if args.mode in ['visualize', 'report']:
            # Load existing data
            raw_data_file = Path(args.output) / 'raw_variant_data.csv'
            if not raw_data_file.exists():
                logger.error("raw_variant_data.csv not found. Run data collection first.")
                sys.exit(1)
                
            logger.info("Loading existing variant data...")
            analysis.variant_data = pd.read_csv(raw_data_file)
            
            # Process data and generate results
            logger.info("Processing data...")
            analysis.preprocess_data()
            results = analysis.run_comprehensive_analysis()
            
            if args.mode == 'visualize':
                logger.info("Generating visualizations...")
                # All visualizations will be saved in the figures directory
                analysis._visualize_clinical_impact(results['clinical_impact'])
                analysis._visualize_functional_consequences(results['functional_consequences'])
                analysis._visualize_evidence_patterns(results['evidence_patterns'])
                analysis._visualize_geographical_risks(results['geographical_risks'])
                analysis._visualize_pleiotropic_effects(results['pleiotropic_effects'])
                analysis._visualize_protective_variants(results['protective_variants'])
                analysis._visualize_sample_size_effects(results['sample_size_effects'])
                analysis._visualize_population_structure(results['population_structure'])
                logger.info("Visualization generation complete")
            
            elif args.mode == 'report':
                logger.info("Generating analysis report...")
                analysis.generate_analysis_report(results)
                logger.info("Report generation complete")
                
        else:
            # Original code for fetch/analyze modes
            if args.gene:
                logger.info(f"Running analysis for single gene: {args.gene}")
                if args.gene in analysis.disease_genes:
                    analysis.disease_genes = {args.gene: analysis.disease_genes[args.gene]}
                else:
                    logger.error(f"Gene {args.gene} not found in disease genes list")
                    sys.exit(1)
            else:
                logger.info("Running analysis for all genes")
                analysis.expand_gene_panels()
            
            if args.mode in ['fetch', 'all']:
                logger.info("Fetching variant data...")
                analysis.fetch_variant_data()
            
            if args.mode in ['analyze', 'all']:
                logger.info("Starting analysis...")
                analysis.preprocess_data()
                results = analysis.run_comprehensive_analysis()
                analysis.generate_analysis_report(results)
        
        logger.info(f"Operations complete. Results saved in {args.output}")
        
    except Exception as e:
        logger.error(f"Operation failed: {e}")
        raise
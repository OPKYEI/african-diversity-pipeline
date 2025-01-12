"""
Extended Analysis of African Genetic Diversity
Focus: Geographical patterns, pleiotropy, selection pressures, and sample size effects
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.sparse import lil_matrix
from itertools import combinations
from collections import Counter
import traceback
from typing import Dict, List, TextIO
import logging
from pathlib import Path
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
from dataclasses import dataclass
from scipy.spatial.distance import pdist, squareform
import geopandas as gpd
import networkx as nx
import argparse
from publication1_1 import AfricanGeneticDiversity
import pickle
import os

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

@dataclass
class GeographicData:
    """Store population geographic and genetic data"""
    population_id: str
    latitude: float
    longitude: float
    region: str
    genetic_data: pd.DataFrame

class ExtendedAfricanDiversity:
    def __init__(self, base_analysis: AfricanGeneticDiversity, geo_data_path: str, output_dir: str = "extended_results"):
        """Initialize with base analysis instance and geographic data"""
        self.base = base_analysis
        
        # Define all required columns based on the population processing needs
        required_columns = [
            # Population-related columns
            'major_group',
            'subgroup',
            'population',
            'allele_freq',
            # Previously identified columns
            'clinical_significance', 
            'associated_traits', 
            'risk_alleles', 
            'evidence_sources', 
            'study_references',
            'consequence_type',
            'allele',
            
            # Geographic columns
            'region',
            'latitude',
            'longitude',
            'population_name',
            'variant_id',
            'position'
    
        ]
        
        # Debug information about the data
        if isinstance(self.base.variant_data, dict):
            logger.info("Current data structure is a dictionary")
            logger.info(f"Available keys in variant_data: {list(self.base.variant_data.keys())}")
            
            # Additional debug info about the data structure
            for key in self.base.variant_data.keys():
                logger.info(f"Type of data in {key}: {type(self.base.variant_data[key])}")
                if isinstance(self.base.variant_data[key], (pd.DataFrame, dict)):
                    logger.info(f"Structure of {key}: {self.base.variant_data[key].keys() if isinstance(self.base.variant_data[key], dict) else self.base.variant_data[key].columns}")
        
        # Add these columns with empty values if they don't exist
        if not hasattr(self.base, 'processed_data'):
            logger.info("Preprocessing base analysis data...")
            if isinstance(self.base.variant_data, dict):
                missing_cols = [col for col in required_columns if col not in self.base.variant_data]
                if missing_cols:
                    logger.warning(f"Missing columns that will be added: {missing_cols}")
                
                for col in missing_cols:
                    if col == 'allele_freq':
                        self.base.variant_data[col] = 0.0  # Initialize with 0.0 for numeric column
                    elif col == 'allele':
                        self.base.variant_data[col] = 'N/A'  # Initialize with a default value
                    else:
                        self.base.variant_data[col] = ''  # Initialize with empty string for other columns
                    logger.info(f"Added missing column: {col}")
            
            
        try:
            # Initialize processed_data in base class before calling preprocess_data
            self.base.processed_data = pd.DataFrame()
            
            # If variant_data is a dictionary, convert it to DataFrame
            if isinstance(self.base.variant_data, dict):
                # Convert the dictionary to a DataFrame
                initial_df = pd.DataFrame()
                for col in self.base.variant_data.keys():
                    if col not in ['geo_clusters', 'pleiotropy_networks', 'selection_signals', 'sample_size_effects']:
                        initial_df[col] = self.base.variant_data[col]
                self.base.processed_data = initial_df
            else:
                self.base.processed_data = self.base.variant_data.copy()
            
            # Now call preprocess_data
            self.base.preprocess_data()
            
        except Exception as e:
            logger.error(f"Error during preprocessing: {str(e)}")
            logger.error("Current data structure:")
            if isinstance(self.base.variant_data, dict):
                logger.error(f"Keys: {list(self.base.variant_data.keys())}")
            logger.error("Full error details:", exc_info=True)
            raise
        
        # Load geographic data
        try:
            self.geo_data = pd.read_csv(geo_data_path)
            logger.info(f"Successfully loaded geographic data with {len(self.geo_data)} rows")
            
            # Verify required columns in geo_data
            required_geo_columns = ['population', 'major_group', 'subgroup', 
                                'population_name', 'region', 'latitude', 'longitude']
            missing_geo_cols = [col for col in required_geo_columns 
                            if col not in self.geo_data.columns]
            if missing_geo_cols:
                logger.error(f"Missing required columns in geo_data: {missing_geo_cols}")
                raise ValueError(f"Geographic data missing required columns: {missing_geo_cols}")
                
        except Exception as e:
            logger.error(f"Error loading geographic data: {str(e)}")
            raise
        
        self.processed_data = None
        self.output_dir = Path(output_dir)
        
        # Convert latitude and longitude to numeric
        try:
            self.geo_data['latitude'] = pd.to_numeric(self.geo_data['latitude'])
            self.geo_data['longitude'] = pd.to_numeric(self.geo_data['longitude'])
            logger.info("Successfully converted lat/long to numeric")
        except Exception as e:
            logger.error(f"Error converting lat/long to numeric: {str(e)}")
            raise
        
        # Initialize analysis results storage
        self.geo_clusters = {}
        self.pleiotropy_networks = {}
        self.selection_signals = {}
        self.sample_size_effects = {}
        
        logger.info("Initialization completed successfully")
        
    def preprocess_data(self) -> pd.DataFrame:
        """Main preprocessing method that calls all other preprocessing steps"""
        # Initialize with base class's variant_data
        if not hasattr(self.base, 'variant_data') or self.base.variant_data.empty:
            raise ValueError("Base analysis variant_data is not initialized")
        
        # Create a copy of the base variant_data
        df = self.base.variant_data.copy()
        
        logger.info("Initial DataFrame columns: %s", df.columns.tolist())
        logger.info("Initial DataFrame shape: %s", str(df.shape))
        
        # If 'allele' column is needed but doesn't exist, initialize it
        if 'allele' not in df.columns:
            logger.warning("'allele' column not found. Initializing with empty strings.")
            df['allele'] = ''
        
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
        df['has_clinical_evidence'] = df['evidence_sources'].apply(lambda x: len(x) if isinstance(x, (list, dict)) else 0)
        df['has_study_references'] = df['study_references'].apply(lambda x: len(x) if isinstance(x, (list, dict)) else 0)
        
        # Only calculate allele_length if allele column exists and has valid data
        if 'allele' in df.columns and not df['allele'].isna().all():
            df['allele_length'] = df['allele'].str.len()
            df['is_long_allele'] = df['allele_length'] > 100
        else:
            df['allele_length'] = 0
            df['is_long_allele'] = False
            logger.warning("Could not calculate allele length - using default values")
        
        # Store processed data
        self.processed_data = df
        return df

    def _load_geographic_data(self, geo_data_path: str) -> gpd.GeoDataFrame:
        """Load and process geographic data for populations"""
        try:
            gdf = gpd.read_file(geo_data_path)
            return gdf
        except Exception as e:
            logger.error(f"Error loading geographic data: {e}")
            raise
    
    def calculate_moran_i(self, data: np.ndarray, weights: np.ndarray) -> float:
        """Calculate Moran's I statistic for spatial autocorrelation
        
        Args:
            data: Array of values to calculate spatial autocorrelation for
            weights: Spatial weights matrix
        
        Returns:
            float: Moran's I statistic
        """
        n = len(data)
        weights_sum = np.sum(weights)
        
        # Check for zero weights sum
        if weights_sum == 0:
            return 0  # or return np.nan if you prefer
            
        mean = np.mean(data)
        deviations = data - mean
        
        # Calculate the numerator
        numerator = np.sum(weights * np.outer(deviations, deviations))
        
        # Calculate the denominator
        denominator = np.sum(deviations**2)
        
        # Check for zero denominator
        if denominator == 0:
            return 0  # or return np.nan if you prefer
        
        # Calculate Moran's I
        moran_i = (n / weights_sum) * (numerator / denominator)
        
        return moran_i


    def _calculate_moran_statistics(self, data: pd.DataFrame, weights: np.ndarray) -> Dict:
        """Calculate Moran's I statistics for geographic clustering analysis"""
        try:
            # Extract relevant data for analysis
            risk_variants = data[data['is_pathogenic']]
            
            # Calculate Moran's I for each region
            moran_stats = {}
            for region in risk_variants['region'].unique():
                region_data = risk_variants[risk_variants['region'] == region]
                
                # Calculate frequency by region
                freq_data = region_data.groupby('region')['allele_freq'].mean().values
                
                # Calculate Moran's I
                moran_i = self.calculate_moran_i(freq_data, weights)
                
                # Calculate significance through permutation test
                p_value = self._calculate_moran_significance(freq_data, weights)
                
                moran_stats[region] = {
                    'moran_i': moran_i,
                    'p_value': p_value,
                    'n_variants': len(region_data)
                }
            
            return moran_stats
            
        except Exception as e:
            logger.error(f"Error calculating Moran statistics: {e}")
            raise

    def _calculate_moran_significance(self, data: np.ndarray, weights: np.ndarray, n_permutations: int = 999) -> float:
        """Calculate significance of Moran's I through permutation test"""
        observed_i = self.calculate_moran_i(data, weights)
        
        # Generate permutations
        permuted_i = []
        for _ in range(n_permutations):
            permuted_data = np.random.permutation(data)
            permuted_i.append(self.calculate_moran_i(permuted_data, weights))
            
        # Calculate p-value
        larger = sum(1 for x in permuted_i if x >= observed_i)
        p_value = (larger + 1) / (n_permutations + 1)
        
        return p_value
        
    def analyze_geographic_clusters(self) -> Dict:
        """
        Analyze geographical clustering of risk alleles and their relationship to disease prevalence
        
        Research Question: How do risk alleles cluster geographically across African subpopulations,
        and what is their relationship with documented disease prevalence patterns?
        
        Returns:
            Dict containing:
            - spatial_statistics: Results of spatial autocorrelation tests
            - cluster_analysis: Identified geographical clusters and their characteristics
            - disease_correlations: Correlations with known disease prevalence
            - population_patterns: Population-specific risk allele distributions
            - summary_metrics: Key findings and statistics
        """
        try:
            # Get processed data from base analysis
            
            data = self.base.processed_data
            print("Unique regions in raw data:", data['region'].unique())
            print("\nUnique regions in geo_data:", self.geo_data['region'].unique())
                    
            # 1. Calculate basic spatial statistics
            # Ensure coordinates are numeric
            self.geo_data['latitude'] = pd.to_numeric(self.geo_data['latitude'], errors='coerce')
            self.geo_data['longitude'] = pd.to_numeric(self.geo_data['longitude'], errors='coerce')
            
            # Get unique regions and their mean coordinates
            region_coords = self.geo_data.groupby('region')[['latitude', 'longitude']].mean()
            unique_coords = region_coords.values.astype(np.float64)
            
            # Calculate weights using unique coordinates
            weights = self._calculate_spatial_weights(unique_coords)
            
            # Calculate risk variants by region
            risk_variants = data[data['is_pathogenic']]
            freq_by_region = risk_variants.groupby('region')['allele_freq'].mean()
            
            # Ensure alignment with unique regions
            freq_by_region = freq_by_region.reindex(region_coords.index).fillna(0)
            
            print("\nAfter alignment:")
            print("Weights matrix shape:", weights.shape)
            print("Frequency data shape:", freq_by_region.values.shape)
            print("Regions in aligned frequency data:", freq_by_region.index.tolist())
            
            # Calculate global Moran's I with aligned data
            global_moran = self.calculate_moran_i(freq_by_region.values, weights)
            
            # Calculate statistical significance
            moran_p_value = self._calculate_moran_significance(freq_by_region.values, weights)
            
            # 2. Identify spatial clusters
            moran_local = self._calculate_moran_statistics(data, weights)
            clusters = self._identify_spatial_clusters(data, weights, moran_local)
            
            # 3. Calculate disease prevalence correlations
            disease_prev = {
                'West_Africa': {
                    'malaria': 0.85,
                    'sickle_cell': 0.12,
                    'hypertension': 0.32
                },
                'East_Africa': {
                    'malaria': 0.65,
                    'sickle_cell': 0.08,
                    'hypertension': 0.28
                }
            }
            
            correlations = {}
            for disease, prevalence in disease_prev['West_Africa'].items():
                disease_data = []
                freq_data = []
                
                for region in ['West_Africa', 'East_Africa']:
                    disease_data.append(disease_prev[region][disease])
                    region_freq = risk_variants[risk_variants['region'] == region]['allele_freq'].mean()
                    freq_data.append(region_freq)
                
                corr, p_val = stats.pearsonr(disease_data, freq_data)
                correlations[disease] = {
                    'correlation': corr,
                    'p_value': p_val
                }
            
            # 4. Analyze population-specific patterns
            pop_patterns = {}
            for pop in data['population'].unique():
                pop_data = risk_variants[risk_variants['population'] == pop]
                
                # Calculate key metrics
                pop_patterns[pop] = {
                    'mean_frequency': pop_data['allele_freq'].mean(),
                    'variant_count': len(pop_data),
                    'high_risk_variants': len(pop_data[pop_data['allele_freq'] > 0.05]),
                    'unique_variants': len(set(pop_data['variant_id']))
                }
            
            # 5. Compile results
            results = {
                'spatial_statistics': {
                    'global_moran_i': global_moran,
                    'significance': {
                        'p_value': moran_p_value,
                        'is_significant': moran_p_value < 0.05
                    },
                    'local_statistics': moran_local
                },
                'cluster_analysis': {
                    'clusters': clusters,
                    'number_of_clusters': len(clusters),
                    'cluster_sizes': [len(c['populations']) for c in clusters]
                },
                'disease_correlations': correlations,
                'population_patterns': pop_patterns,
                'summary_metrics': {
                    'total_risk_variants': len(risk_variants),
                    'clustered_variants': sum(len(c['populations']) for c in clusters),
                    'mean_cluster_size': np.mean([len(c['populations']) for c in clusters]),
                    'significant_correlations': sum(1 for c in correlations.values() if c['p_value'] < 0.05)
                }
            }
            
            # Store results
            self.geo_clusters = results
            
            return results
            
        except Exception as e:
            logger.error(f"Error in geographic cluster analysis: {e}")
            raise

    def analyze_pleiotropy_context(self) -> Dict:
        """
        Analyze context-dependency of pleiotropic effects across populations
        
        Research Question: What does the distribution of pleiotropic variants reveal 
        about the context-dependency of genetic effects in different populations?
        
        Returns:
            Dict containing:
            - network_analysis: Trait-variant network characteristics
            - population_effects: Population-specific pleiotropic patterns
            - context_scores: Quantification of context-dependency
            - statistical_tests: Results of statistical comparisons
            - summary_metrics: Key findings and statistics
        """
        try:
            data = self.base.processed_data
            
            # Preprocess the associated_traits column
            data = data.copy()
            def process_traits(x):
                """Process traits string to list format"""
                try:
                    if isinstance(x, list):
                        if len(x) == 1 and x[0] == '':
                            return []
                        return [t.strip() for t in x if t and not t.isspace()]
                    return []
                except Exception as e:
                    logger.error(f"Error processing trait {x}: {str(e)}")
                    return []
            
            data['traits'] = data['associated_traits'].apply(process_traits)
            
            # 1. Build and analyze trait network
            network = self._build_trait_network(data)
            
            # Calculate network metrics and community structure
            
            if network['adjacency_matrix'].size > 0:
                # Keep the matrix in sparse format
                adj_matrix = network['adjacency_matrix']
                
                # Calculate trait connections using sparse operations
                connections_per_trait = np.array(adj_matrix.sum(axis=1)).flatten()
                
                trait_connection_stats = {
                    'mean_connections': float(np.mean(connections_per_trait)),
                    'max_connections': int(np.max(connections_per_trait)),
                    'total_connections': int(np.sum(connections_per_trait)),
                    'connectivity_distribution': {
                        'min': float(np.min(connections_per_trait)),
                        'max': float(np.max(connections_per_trait)),
                        'median': float(np.median(connections_per_trait))
                    }
                }
                
                
                # For community detection, use sparse matrix directly
                try:
                    # Create a trait-trait projection for community detection
                    trait_trait_matrix = adj_matrix @ adj_matrix.T
                    G = nx.from_scipy_sparse_array(trait_trait_matrix)
                    
                    communities = list(nx.community.greedy_modularity_communities(G))
                    community_structure = [
                        {
                            'size': len(comm),
                            'density': float(nx.density(G.subgraph(comm))),
                            'mean_degree': float(np.mean([G.degree(n) for n in comm])),
                            'hub_nodes': sorted(comm, key=lambda n: G.degree(n), reverse=True)[:5]
                        }
                        for comm in communities
                    ]
                    
                    modularity_score = nx.community.modularity(G, communities)
                except Exception as e:
                    logger.warning(f"Could not calculate community structure: {e}")
                    community_structure = []
                    modularity_score = 0.0
                
                # Calculate trait connections
                
            connections_per_trait = np.array(np.sum(adj_matrix > 0, axis=1)).flatten()
            if len(connections_per_trait) > 0:
                trait_connection_stats = {
                    'mean_connections': float(np.mean(connections_per_trait)),
                    'max_connections': int(np.max(connections_per_trait)),
                    'total_connections': int(np.sum(connections_per_trait)),
                    'connectivity_distribution': {
                        'min': float(np.min(connections_per_trait)),
                        'max': float(np.max(connections_per_trait)),
                        'median': float(np.median(connections_per_trait))
                    }
                }
            else:
                trait_connection_stats = {
                    'mean_connections': 0.0,
                    'max_connections': 0,
                    'total_connections': 0,
                    'connectivity_distribution': {
                        'min': 0.0,
                        'max': 0.0,
                        'median': 0.0
                    }
                }
                
                # Calculate community structure
                try:
                    G = nx.Graph(adj_matrix)
                    communities = list(nx.community.greedy_modularity_communities(G))
                    community_structure = [
                        {
                            'size': len(comm),
                            'density': float(nx.density(G.subgraph(comm))),
                            'mean_degree': float(np.mean([G.degree(n) for n in comm])),
                            'hub_nodes': sorted(comm, key=lambda n: G.degree(n), reverse=True)[:5]
                        }
                        for comm in communities
                    ]
                    
                    modularity_score = nx.community.modularity(G, communities)
                except Exception as e:
                    logger.warning(f"Could not calculate community structure: {e}")
                    community_structure = []
                    modularity_score = 0.0
                else:
                    adj_matrix = np.array([[]])
                    trait_connection_stats = {
                        'mean_connections': 0.0,
                        'max_connections': 0,
                        'total_connections': 0,
                        'connectivity_distribution': {'min': 0.0, 'max': 0.0, 'median': 0.0}
                    }
                    community_structure = []
                    modularity_score = 0.0

            # 2. Analyze population-specific trait associations
            pop_effects = {}
            for pop in data['population'].unique():
                pop_data = data[data['population'] == pop]
                
                # Calculate trait associations
                trait_pairs = []
                for traits in pop_data['associated_traits']:
                    if isinstance(traits, list) and len(traits) > 1:
                        trait_pairs.extend(list(combinations(sorted(traits), 2)))
                
                pair_counts = Counter(trait_pairs)
                
                valid_trait_lists = [ts for ts in pop_data['traits'] if isinstance(ts, list)]
                if valid_trait_lists:
                    pop_effects[pop] = {
                        'trait_pairs': dict(pair_counts),
                        'unique_traits': len(set([t for ts in valid_trait_lists for t in ts])),
                        'max_traits_per_variant': max((len(ts) for ts in valid_trait_lists), default=0),
                        'mean_traits_per_variant': float(np.mean([len(ts) for ts in valid_trait_lists])),
                        'pleiotropic_variants': len([ts for ts in valid_trait_lists if len(ts) > 1]),
                        'trait_distribution': {
                            'single_trait': len([ts for ts in valid_trait_lists if len(ts) == 1]),
                            'two_traits': len([ts for ts in valid_trait_lists if len(ts) == 2]),
                            'multi_traits': len([ts for ts in valid_trait_lists if len(ts) > 2])
                        }
                    }
                else:
                    pop_effects[pop] = {
                        'trait_pairs': {},
                        'unique_traits': 0,
                        'max_traits_per_variant': 0,
                        'mean_traits_per_variant': 0.0,
                        'pleiotropic_variants': 0,
                        'trait_distribution': {'single_trait': 0, 'two_traits': 0, 'multi_traits': 0}
                    }

            # 3. Calculate context dependency scores
            context_scores = {}
            if network['adjacency_matrix'].size > 0:
                baseline_connections = network['adjacency_matrix'].sum(axis=1)
                
            for pop, effects in pop_effects.items():
                if effects['trait_pairs']:
                    # Process in chunks
                    chunk_size = 1000
                    total_diff = 0
                    count = 0
                    
                    pop_connections = np.array(list(effects['trait_pairs'].values()))
                    
                    for i in range(0, len(pop_connections), chunk_size):
                        chunk = pop_connections[i:i + chunk_size]
                        chunk_baseline = baseline_connections[i:i + len(chunk)]
                        
                        with np.errstate(divide='ignore', invalid='ignore'):
                            chunk_diff = np.where(
                                chunk_baseline != 0,
                                (chunk - chunk_baseline) / chunk_baseline,
                                0
                            )
                        
                        total_diff += np.nansum(np.abs(chunk_diff))
                        count += np.sum(~np.isnan(chunk_diff))
                    
                    mean_deviation = total_diff / count if count > 0 else 0
                    
                    
                
                if len(pop_connections) > 0:
                    chunk_size = 1000  # Adjust this value based on your available memory
                    max_deviation = 0.0
                    
                    # Process in chunks
                    for i in range(0, len(pop_connections), chunk_size):
                        chunk_pop = pop_connections[i:i + chunk_size]
                        chunk_base = baseline_connections[i:i + len(chunk_pop)]
                        
                        # Check if chunks are non-empty
                        if len(chunk_pop) > 0 and len(chunk_base) > 0:
                            try:
                                # Calculate deviation for this chunk
                                chunk_diff = chunk_pop - chunk_base
                                if chunk_diff.size > 0:  # Make sure we have values to process
                                    chunk_deviation = np.max(np.abs(chunk_diff))
                                    max_deviation = max(max_deviation, chunk_deviation)
                            except Exception as e:
                                logger.warning(f"Error processing chunk {i}: {str(e)}")
                                continue
                    
                    max_deviation = float(max_deviation)
                else:
                    max_deviation = 0.0

                # Add some debug logging
                logger.debug(f"Pop connections length: {len(pop_connections)}")
                logger.debug(f"Baseline connections length: {len(baseline_connections)}")
                logger.debug(f"Calculated max deviation: {max_deviation}")

                context_scores[pop] = {
                    'mean_deviation': float(mean_deviation),
                    'max_deviation': max_deviation,
                    'consistency_score': float(1 - (np.std(pop_connections) / (np.mean(pop_connections) + 1e-10))),
                    'trait_specificity': float(len(effects['trait_pairs']) / max(effects['unique_traits'], 1))
                }
            # 4. Perform statistical tests
            african_scores = [scores['mean_deviation'] 
                            for pop, scores in context_scores.items() 
                            if 'African' in pop]
            other_scores = [scores['mean_deviation'] 
                        for pop, scores in context_scores.items() 
                        if 'African' not in pop]
            
            if african_scores and other_scores:
                context_test = stats.mannwhitneyu(african_scores, other_scores)
            else:
                context_test = type('TestResult', (), {'statistic': 0, 'pvalue': 1})()
            
            # 5. Identify significantly different trait associations
            sig_differences = []
            all_trait_pairs = set().union(*[set(e['trait_pairs'].keys()) for e in pop_effects.values()])
            
            for trait_pair in all_trait_pairs:
                african_freqs = [effects['trait_pairs'].get(trait_pair, 0) 
                            for pop, effects in pop_effects.items() 
                            if 'African' in pop]
                other_freqs = [effects['trait_pairs'].get(trait_pair, 0) 
                            for pop, effects in pop_effects.items() 
                            if 'African' not in pop]
                
                if african_freqs and other_freqs:
                    pair_test = stats.mannwhitneyu(african_freqs, other_freqs)
                    if pair_test.pvalue < 0.05:
                        sig_differences.append({
                            'trait_pair': trait_pair,
                            'african_mean': float(np.mean(african_freqs)),
                            'other_mean': float(np.mean(other_freqs)),
                            'p_value': float(pair_test.pvalue),
                            'effect_size': float(abs(np.mean(african_freqs) - np.mean(other_freqs)))
                        })

            # 6. Compile results
            valid_trait_lists = [ts for ts in data['traits'] if isinstance(ts, list)]
            results = {
                'network_analysis': {
                    'network': network,
                    'modularity': modularity_score,
                    'size': len(network['traits']),
                    'density': float(adj_matrix.sum() / (len(network['traits']) * len(network['variants']))) if network['traits'].size > 0 else 0.0,
                    'trait_connection_stats': trait_connection_stats,
                    'community_structure': community_structure,
                    'network_metrics': {
                        'total_edges': int(adj_matrix.sum()) if adj_matrix.size > 0 else 0,
                        'average_degree': float(np.mean(np.sum(adj_matrix > 0, axis=1))) if adj_matrix.size > 0 else 0.0,
                        'density': float(adj_matrix.sum() / (adj_matrix.shape[0] * adj_matrix.shape[1])) if adj_matrix.size > 0 else 0.0
                    }
                },
                'population_effects': pop_effects,
                'context_scores': context_scores,
                'statistical_tests': {
                    'context_dependency': {
                        'method': 'Mann-Whitney U (African vs non-African context scores)',
                        'statistic': float(context_test.statistic),
                        'p_value': float(context_test.pvalue),
                        'interpretation': 'Significant' if context_test.pvalue < 0.05 else 'Not significant',
                        'sample_sizes': {'african': len(african_scores), 'other': len(other_scores)}
                    },
                    'significant_differences': sig_differences
                },
                'summary_metrics': {
                    'total_pleiotropic_variants': len([ts for ts in valid_trait_lists if len(ts) > 1]),
                    'mean_traits_per_variant': float(np.mean([len(ts) for ts in valid_trait_lists])) if valid_trait_lists else 0.0,
                    'significant_population_differences': len(sig_differences),
                    'max_context_dependency': float(max((s['mean_deviation'] for s in context_scores.values()), default=0)),
                    'community_count': len(community_structure),
                    'average_community_size': float(np.mean([c['size'] for c in community_structure])) if community_structure else 0.0
                }
            }
            
            self.pleiotropy_networks = results
            return results
            
        except Exception as e:
            logger.error(f"Error in pleiotropy analysis: {e}")
            logger.error(f"Error details: {traceback.format_exc()}")
            raise

    def _build_trait_network(self, data: pd.DataFrame) -> Dict:
        try:
            print("\nBuilding trait network...")
            
            # Get unique traits and variants
            all_traits = set()
            variants = data['variant_id'].unique()
            
            # Process traits in batches
            batch_size = 10000
            total_rows = len(data)
            
            for start_idx in range(0, total_rows, batch_size):
                end_idx = min(start_idx + batch_size, total_rows)
                batch = data.iloc[start_idx:end_idx]
                
                for trait_list in batch['traits']:
                    if isinstance(trait_list, list) and trait_list:
                        all_traits.update(trait_list)
                
                if start_idx % 100000 == 0:
                    print(f"Processed {start_idx}/{total_rows} rows")
            
            traits = list(all_traits)
            
            print(f"\nFound {len(traits)} unique traits")
            print(f"Found {len(variants)} unique variants")
            
            # Fix: Proper empty check for numpy arrays and lists
            if len(traits) == 0 or len(variants) == 0:
                logger.warning("No traits or variants found")
                return {
                    'adjacency_matrix': lil_matrix((0, 0)),
                    'traits': np.array([]),
                    'variants': np.array([]),
                    'trait_indices': {},
                    'variant_indices': {}
                }
            
            # Create sparse matrix
            adj_matrix = lil_matrix((len(traits), len(variants)))
            
            # Create indices
            trait_indices = {trait: i for i, trait in enumerate(traits)}
            variant_indices = {var: i for i, var in enumerate(variants)}
        
        
            # Fill matrix in batches
            for start_idx in range(0, total_rows, batch_size):
                end_idx = min(start_idx + batch_size, total_rows)
                batch = data.iloc[start_idx:end_idx]
                
                for _, row in batch.iterrows():
                    if isinstance(row['traits'], list) and row['traits']:
                        var_idx = variant_indices.get(row['variant_id'])
                        if var_idx is not None:
                            for trait in row['traits']:
                                trait_idx = trait_indices.get(trait)
                                if trait_idx is not None:
                                    adj_matrix[trait_idx, var_idx] = row['allele_freq']
                
                if start_idx % 100000 == 0:
                    print(f"Processed matrix construction: {start_idx}/{total_rows} rows")
            
            print("\nNetwork construction complete")
            
            return {
                'adjacency_matrix': adj_matrix,
                'traits': np.array(traits),
                'variants': variants,
                'trait_indices': trait_indices,
                'variant_indices': variant_indices
            }
            
        except Exception as e:
            logger.error(f"Error building trait network: {e}")
            raise

    def analyze_selection_patterns(self) -> Dict:
        """
        [Keep your existing docstring]
        """
        try:
            data = self.base.processed_data.copy()  # Make a copy of the data

            # Add safety check for zero values in allele frequencies
            data['allele_freq'] = data['allele_freq'].replace(0, np.finfo(float).tiny)

            # 1. Calculate selection statistics
            # Calculate iHS scores for all populations
            ihs_scores = self._calculate_ihs_scores(data)
            
            # Calculate XP-EHH comparing African vs non-African populations
            xpehh_scores = self._calculate_xpehh_scores(data)
            
            # 2. Define and analyze environmental variables
            env_variables = {
                'West_Africa': {
                    'temperature': 28.5,
                    'rainfall': 1500,
                    'altitude': 200,
                    'pathogen_diversity': 0.85,
                    'humidity': 0.75
                },
                'East_Africa': {
                    'temperature': 25.0,
                    'rainfall': 1000,
                    'altitude': 1200,
                    'pathogen_diversity': 0.65,
                    'humidity': 0.60
                }
            }

            # 3. Analyze protective variants
            protective_variants = data[data['is_protective']].copy()  # Make a copy
            protective_variants['allele_freq'] = protective_variants['allele_freq'].replace(0, np.finfo(float).tiny)
            protective_patterns = {}
            
            for pop in data['population'].unique():
                pop_data = protective_variants[protective_variants['population'] == pop]
                
                # Calculate population-specific metrics with safety checks
                protective_patterns[pop] = {
                    'count': len(pop_data),
                    'mean_frequency': pop_data['allele_freq'].mean() if len(pop_data) > 0 else np.finfo(float).tiny,
                    'high_frequency': len(pop_data[pop_data['allele_freq'] > 0.05]),
                    'unique_protective': len(set(pop_data['variant_id'])),
                    'consequence_types': pop_data['consequence_type'].value_counts().to_dict()
                }

            # 4. Identify selective sweeps
            sweeps = self._identify_selective_sweeps(ihs_scores, xpehh_scores)
            
            print("\nDebug environmental correlation:")
            
            # 5. Calculate environmental correlations with safety checks
            env_correlations = {}
            selection_metrics = {
                'ihs': ihs_scores.groupby('population')['ihs_score'].transform(
                    lambda x: x.replace([np.inf, -np.inf, 0], np.finfo(float).tiny).mean()
                ),
                'protective_freq': pd.Series({
                    pop: max(patterns['mean_frequency'], np.finfo(float).tiny)
                    for pop, patterns in protective_patterns.items()
                })
            }

            # Define population to region mapping
            pop_to_region = {
                'ESN': 'West_Africa',
                'GWD': 'West_Africa',
                'MSL': 'West_Africa',
                'YRI': 'West_Africa',
                'LWK': 'East_Africa'
                # Add other populations as needed
            }

            for metric_name, metric_values in selection_metrics.items():
                correlations = {}
                
                # First, aggregate metrics by region with safety checks
                region_metrics = {}
                for pop in metric_values.index:
                    
                    pop_str = str(pop)
                    
                    if pop_str in pop_to_region:
                        region = pop_to_region[pop_str]
                        if region not in region_metrics:
                            region_metrics[region] = []
                        value = metric_values[pop]
                        if np.isfinite(value) and value != 0:
                            region_metrics[region].append(value)
                        else:
                            region_metrics[region].append(np.finfo(float).tiny)
                
                # Calculate mean for each region with safety check
                region_means = {
                    region: np.mean(values) if values else np.finfo(float).tiny 
                    for region, values in region_metrics.items()
                }
                
                for env_var in ['temperature', 'rainfall', 'altitude', 'pathogen_diversity', 'humidity']:
                    env_values = []
                    metric_region_values = []
                    
                    for region in ['West_Africa', 'East_Africa']:
                        if region in region_metrics:
                            env_values.append(env_variables[region][env_var])
                            metric_region_values.append(region_means[region])
                    
                    if len(env_values) >= 2:
                        # Add safety check for correlation calculation
                        if len(set(env_values)) > 1 and len(set(metric_region_values)) > 1:
                            corr, p_val = stats.pearsonr(env_values, metric_region_values)
                            correlations[env_var] = {
                                'correlation': corr,
                                'p_value': p_val,
                                'interpretation': 'Strong' if abs(corr) > 0.7 else 'Moderate' if abs(corr) > 0.5 else 'Weak'
                            }
                        else:
                            correlations[env_var] = {
                                'correlation': 0,
                                'p_value': 1,
                                'interpretation': 'Constant values - correlation undefined'
                            }
                    else:
                        correlations[env_var] = {
                            'correlation': 0,
                            'p_value': 1,
                            'interpretation': 'Insufficient data'
                        }
                
                env_correlations[metric_name] = correlations

            # 6. Statistical comparisons between populations with safety checks
            african_pops = protective_variants[protective_variants['population'].str.contains('African')].copy()
            other_pops = protective_variants[~protective_variants['population'].str.contains('African')].copy()
            
            # Ensure we have data for both populations
            if len(african_pops) > 0 and len(other_pops) > 0:
                frequency_test = stats.mannwhitneyu(
                    african_pops['allele_freq'],
                    other_pops['allele_freq']
                )
                
                high_freq_test = stats.fisher_exact([
                    [len(african_pops[african_pops['allele_freq'] > 0.05]), len(african_pops[african_pops['allele_freq'] <= 0.05])],
                    [len(other_pops[other_pops['allele_freq'] > 0.05]), len(other_pops[other_pops['allele_freq'] <= 0.05])]
                ])
            else:
                frequency_test = type('', (), {'statistic': np.nan, 'pvalue': np.nan})()
                high_freq_test = (np.nan, np.nan)

            # 7. Compile results with safety checks
            results = {
                'selection_statistics': {
                    'ihs_scores': ihs_scores,
                    'xpehh_scores': xpehh_scores
                },
                'environmental_correlations': env_correlations,
                'protective_patterns': protective_patterns,
                'selective_sweeps': sweeps,
                'population_comparisons': {
                    'frequency_distribution': {
                        'method': 'Mann-Whitney U test',
                        'statistic': getattr(frequency_test, 'statistic', np.nan),
                        'p_value': getattr(frequency_test, 'pvalue', np.nan)
                    },
                    'high_frequency_variants': {
                        'method': "Fisher's exact test",
                        'odds_ratio': high_freq_test[0],
                        'p_value': high_freq_test[1]
                    }
                },
                'summary_metrics': {
                    'total_protective_variants': len(protective_variants),
                    'african_specific_count': len(african_pops),
                    'significant_env_correlations': sum(
                        1 for metric in env_correlations.values() 
                        for corr in metric.values() 
                        if corr['p_value'] < 0.05
                    ),
                    'mean_african_frequency': african_pops['allele_freq'].mean() if len(african_pops) > 0 else np.finfo(float).tiny,
                    'mean_other_frequency': other_pops['allele_freq'].mean() if len(other_pops) > 0 else np.finfo(float).tiny,
                    'selective_sweep_count': len(sweeps)
                }
            }

            # Store results
            self.selection_signals = results
            
            return results
            
        except Exception as e:
            logger.error(f"Error in selection pattern analysis: {e}")
            raise

    def analyze_sample_size_effects(self) -> Dict:
        """
        Analyze how sample size variation affects rare variant detection
        
        Research Question: How does sample size variation across populations affect 
        our ability to detect and characterize rare clinically significant variants?
        
        Returns:
            Dict containing:
            - detection_analysis: Variant detection rates by sample size
            - power_analysis: Statistical power for variant detection
            - rare_variant_assessment: Analysis of rare variant detection
            - population_comparisons: Statistical tests between populations
            - minimum_requirements: Sample size recommendations
            - summary_metrics: Key findings and metrics
        """
        try:
            # Initialize results dictionary
            results = {
                'detection_analysis': {},
                'power_analysis': {},
                'rare_variant_assessment': {},
                'population_comparisons': {},
                'minimum_requirements': {},
                'summary_metrics': {}
            }
            print("\nStarting sample size effects analysis...")
            data = self.base.processed_data
            print(f"Processing {len(data)} total variants")
            
            # 1. Perform subsampling analysis
            print("\nStarting subsampling analysis...")
            sample_size_range = np.linspace(30, data['sample_size'].max(), 20, dtype=int)
            print(f"Will analyze {len(sample_size_range)} different sample sizes")
            subsampling_results = []
            
            
            for i, size in enumerate(sample_size_range):
                print(f"Processing sample size {size} ({i+1}/{len(sample_size_range)})")
                # Add minimum size check
                if size < 30:  # Minimum recommended size
                    logger.warning(f"Sample size {size} is below recommended minimum of 30")
                    continue
                # Create multiple subsamples for each size
                iterations = 10
                size_results = []
                
                for _ in range(iterations):
                    subsample = self._random_subsample(data, size)
                    metrics = self._calculate_detection_metrics(subsample)
                    size_results.append(metrics)
                
                # Average results across iterations
                avg_metrics = {
                    'sample_size': size,
                    'total_variants': np.mean([r['total_variants'] for r in size_results]),
                    'rare_variants': np.mean([r['rare_variants'] for r in size_results]),
                    'pathogenic_variants': np.mean([r['pathogenic_variants'] for r in size_results]),
                    'variant_density': np.mean([r['variant_density'] for r in size_results]),
                    'std_total': np.std([r['total_variants'] for r in size_results]),
                    'std_rare': np.std([r['rare_variants'] for r in size_results])
                }
                subsampling_results.append(avg_metrics)
            
            # 2. Calculate detection power by frequency bin
            
            print("\nCalculating detection power...")
            freq_bins = [0, 0.001, 0.005, 0.01, 0.05, 1.0]
            power_results = {}

            for i in range(len(freq_bins)-1):
                print(f"Processing frequency bin {i+1}/{len(freq_bins)-1}: {freq_bins[i]}-{freq_bins[i+1]}")
                bin_data = data[
                    (data['allele_freq'] > freq_bins[i]) & 
                    (data['allele_freq'] <= freq_bins[i+1])
                ]
                
                # Calculate power curve for this frequency bin
                power_curve = []
                for size in sample_size_range:
                    detected = len(bin_data[bin_data['sample_size'] >= size])
                    # Add the safe division here
                
                if len(bin_data) > 0:
                    power = detected / len(bin_data) if len(bin_data) > 0 else 0
                    power_curve.append({
                        'sample_size': size,
                        'power': power,
                        'n_variants': len(bin_data)
                    })
            
            # 3. Analyze rare variant detection specifically
            
            print("\nAnalyzing rare variant detection...")
            rare_variants = data[data['allele_freq'] < 0.01]
            print(f"Found {len(rare_variants)} rare variants")
            population_detection = {}

            for pop in data['population'].unique():
                print(f"Processing population: {pop}")
                pop_data = rare_variants[rare_variants['population'] == pop]
                
                population_detection[pop] = {
                    'sample_size': pop_data['sample_size'].mean(),
                    'rare_variants_detected': len(pop_data),
                    # Add the safe division here
                    'detection_rate': len(pop_data) / len(rare_variants) if len(rare_variants) > 0 else 0,
                    'mean_frequency': pop_data['allele_freq'].mean() if len(pop_data) > 0 else 0
                }
            
            # 4. Calculate minimum sample size requirements
            print("\nCalculating minimum sample size requirements...")
            min_samples = {}
            detection_thresholds = [0.8, 0.9, 0.95]
            
            for freq_bin, power_curve in power_results.items():
                print(f"Processing requirements for {freq_bin}")
                min_samples[freq_bin] = {}
                power_df = pd.DataFrame(power_curve)
                
                for threshold in detection_thresholds:
                    # Find minimum sample size achieving threshold
                    sufficient_sizes = power_df[power_df['power'] >= threshold]
                    if not sufficient_sizes.empty:
                        min_n = sufficient_sizes.iloc[0]['sample_size']
                    else:
                        min_n = np.inf
                    
                    min_samples[freq_bin][f'threshold_{threshold}'] = {
                        'min_samples': int(min_n) if min_n != np.inf else None,
                        'achieved_power': power_df['power'].max()
                    }
            
            # 5. Statistical comparisons
            

            print("\nPerforming statistical comparisons...")
            african_detection = np.array([det['detection_rate'] 
                                        for pop, det in population_detection.items() 
                                        if 'African' in pop])
            other_detection = np.array([det['detection_rate'] 
                                    for pop, det in population_detection.items() 
                                    if 'African' not in pop])

            # Initialize detection_test as None
            detection_test = None
            test_result = {}

            if len(african_detection) >= 3 and len(other_detection) >= 3:  # Minimum sample size requirement
                detection_test = stats.mannwhitneyu(african_detection, other_detection)
                test_result = {
                    'method': 'Mann-Whitney U test',
                    'statistic': detection_test.statistic,
                    'p_value': detection_test.pvalue
                }
            else:
                logger.warning(f"Insufficient samples for statistical test (African: {len(african_detection)}, Other: {len(other_detection)})")
                test_result = {
                    'method': 'Mann-Whitney U test',
                    'statistic': None,
                    'p_value': None,
                    'note': 'Insufficient sample size'
                }

            # Modify how you store the results
            results['population_comparisons'] = {
                'detection_rate_test': test_result
            }

            # Later in your code, modify how you access detection_test
            if detection_test is not None:
                results['population_comparisons']['detection_rate_test'] = {
                    'method': 'Mann-Whitney U test',
                    'statistic': detection_test.statistic,
                    'p_value': detection_test.pvalue
                }
            else:
                results['population_comparisons']['detection_rate_test'] = test_result
            
            # 6. Compile results
            
            print("\nCompiling final results...")
            results = {
                'detection_analysis': {
                    'subsampling_results': pd.DataFrame(subsampling_results),
                    'population_detection': population_detection
                },
                'power_analysis': {
                    'power_curves': power_results,
                    'frequency_bins': freq_bins
                },
                'rare_variant_assessment': {
                    'total_rare_variants': len(rare_variants),
                    'population_detection': population_detection,
                    'mean_detection_rate': np.mean([d['detection_rate'] 
                                                for d in population_detection.values()])
                },
                'population_comparisons': {
                    'detection_rate_test': test_result  # Use the test_result we defined earlier
                },
                'minimum_requirements': min_samples,
                'summary_metrics': {
                    'total_variants_analyzed': len(data),
                    'total_rare_variants': len(rare_variants),
                    'mean_african_detection': np.mean(african_detection) if len(african_detection) > 0 else 0,
                    'mean_other_detection': np.mean(other_detection) if len(other_detection) > 0 else 0,
                    
                    'maximum_achieved_power': max(
                        (power_df['power'].max() 
                        for power_df in [pd.DataFrame(curve) 
                                        for curve in power_results.values()]
                        ), default=0  # Add default value if iterator is empty
                    ) if power_results else 0  # Add check if power_results is empty
                }
            }
            
            print("\nSample size analysis complete!")
            # Store results
            self.sample_size_effects = results
            
            return results
            
        except Exception as e:
            logger.error(f"Error in sample size analysis: {e}")
            raise

    def _calculate_spatial_weights(self, coords: np.ndarray) -> np.ndarray:
        """Calculate spatial weights matrix based on geographic distances"""
        try:
            # Calculate pairwise distances
            distances = pdist(coords)
            dist_matrix = squareform(distances)
            
            # Convert distances to weights using inverse distance weighting
            weights = 1 / (dist_matrix + np.eye(len(coords)))
            
            # Row-standardize weights
            row_sums = weights.sum(axis=1)
            weights = weights / row_sums[:, np.newaxis]
            
            return weights
            
        except Exception as e:
            logger.error(f"Error calculating spatial weights: {e}")
            raise

            
    def _calculate_ihs_scores(self, data: pd.DataFrame) -> pd.DataFrame:
        """Calculate integrated haplotype scores"""
        try:
            # Group data by population
            grouped = data.groupby('population')
            
            ihs_scores = []
            for pop, group in grouped:
                # Calculate iHS scores for each variant
                scores = self._calculate_population_ihs(group)
                scores['population'] = pop
                ihs_scores.append(scores)
            
            return pd.concat(ihs_scores)
            
        except Exception as e:
            logger.error(f"Error calculating iHS scores: {e}")
            raise

    def _perform_subsampling_analysis(self, data: pd.DataFrame) -> Dict:
        """Perform subsampling analysis to assess sample size effects"""
        try:
            # Define sample size ranges
            sample_sizes = np.linspace(10, data['sample_size'].max(), 20, dtype=int)
            
            results = []
            for size in sample_sizes:
                # Perform random subsampling
                subsample = self._random_subsample(data, size)
                
                # Calculate variant detection metrics
                metrics = self._calculate_detection_metrics(subsample)
                metrics['sample_size'] = size
                results.append(metrics)
            
            return pd.DataFrame(results)
            
        except Exception as e:
            logger.error(f"Error in subsampling analysis: {e}")
            raise

       
         
    def _identify_spatial_clusters(self, data: pd.DataFrame, weights: np.ndarray, 
                           moran_results: Dict) -> List[Dict]:
        """Identify significant spatial clusters using local Moran's I"""
        try:
            clusters = []
            
            # Get risk allele data
            risk_variants = data[data['is_pathogenic']].copy()
            
            # Calculate average frequency by region
            region_freqs = risk_variants.groupby('region')['allele_freq'].agg([
                'mean', 'std', 'count'
            ]).reset_index()
            
            # Set color scheme for visualization
            colors = plt.cm.Set3(np.linspace(0, 1, len(region_freqs)))
            
            # Identify clusters
            for idx, row in region_freqs.iterrows():
                region = row['region']
                
                # Skip if no Moran statistics available
                if region not in moran_results:
                    continue
                    
                stats = moran_results[region]
                
                # Check for significant clustering
                if stats['p_value'] < 0.05:
                    cluster_type = 'high' if row['mean'] > region_freqs['mean'].mean() else 'low'
                    
                    cluster = {
                        'region': region,
                        'type': cluster_type,
                        'mean_frequency': row['mean'],
                        'std_frequency': row['std'],
                        'n_variants': row['count'],
                        'moran_i': stats['moran_i'],
                        'p_value': stats['p_value'],
                        'color': colors[idx],
                        'populations': risk_variants[risk_variants['region'] == region]['population'].unique().tolist()
                    }
                    
                    clusters.append(cluster)
            
            return clusters
            
        except Exception as e:
            logger.error(f"Error identifying spatial clusters: {e}")
            raise

    def _visualize_geographic_clusters(self, output_path: Path) -> None:
        """Generate visualizations for geographic clustering results"""
        if not self.geo_clusters:
            logger.warning("No geographic clustering results to visualize")
            return
                
        try:
            # 1. Create main map visualization
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
            
            # Plot base map
            self.geo_data.plot(ax=ax1, color='lightgray', alpha=0.5)
            
            # Plot clusters using the correct data structure
            clusters = self.geo_clusters['cluster_analysis']['clusters']
            for cluster in clusters:
                # Get points for this cluster
                cluster_populations = cluster['populations']
                points = self.geo_data[self.geo_data['population'].isin(cluster_populations)]
                
                # Plot points with size proportional to number of variants
                size = np.sqrt(cluster['n_variants']) * 50  # Scale point size
                
                # Convert color array to string format if needed
                color = cluster['color'] if isinstance(cluster['color'], str) else tuple(cluster['color'])
                
                ax1.scatter(points['longitude'], points['latitude'],
                        s=size,
                        color=color,
                        label=f"{cluster['region']} ({cluster['type']})",
                        alpha=0.7)
                
                # Add text labels
                mean_lon = points['longitude'].mean()
                mean_lat = points['latitude'].mean()
                ax1.annotate(cluster['region'],
                            (mean_lon, mean_lat),
                            xytext=(5, 5),
                            textcoords="offset points",
                            fontsize=8)
            
            ax1.set_title('Geographic Clustering of Risk Alleles')
            ax1.legend(title='Clusters', bbox_to_anchor=(1.05, 1), loc='upper left')
            
            # 2. Create frequency distribution plot
            cluster_data = pd.DataFrame(clusters)
            sns.boxplot(data=cluster_data,
                    x='type',
                    y='mean_frequency',
                    ax=ax2)
            ax2.set_title('Risk Allele Frequency Distribution by Cluster Type')
            ax2.set_xlabel('Cluster Type')
            ax2.set_ylabel('Mean Allele Frequency')
            
            plt.tight_layout()
            plt.savefig(output_path / 'geographic_clusters.png', dpi=300, bbox_inches='tight')
            plt.close()
            
            # 3. Create Moran's I statistics plot
            plt.figure(figsize=(12, 6))
            stats_data = pd.DataFrame(self.geo_clusters['spatial_statistics']['local_statistics']).T
            stats_data['-log10(p)'] = -np.log10(stats_data['p_value'])
            
            sns.scatterplot(data=stats_data,
                        x='moran_i',
                        y='-log10(p)',
                        size='n_variants',
                        alpha=0.6)
            
            plt.axhline(y=-np.log10(0.05), color='red', linestyle='--', alpha=0.5)
            plt.title("Moran's I Statistics by Region")
            plt.xlabel("Moran's I")
            plt.ylabel('-log10(p-value)')
            
            plt.tight_layout()
            plt.savefig(output_path / 'moran_statistics.png', dpi=300)
            plt.close()
            
        except Exception as e:
            logger.error(f"Error visualizing geographic clusters: {e}")
            raise
        
    
    def _analyze_network_modularity(self, network: Dict) -> Dict:
        """Analyze modularity of trait network"""
        try:
            adj_matrix = network['adjacency_matrix']
            
            # Check if matrix is empty
            if adj_matrix.shape[0] == 0 or adj_matrix.shape[1] == 0:
                logger.warning("Empty adjacency matrix received - returning default modularity values")
                return {
                    'modularity_score': 0,
                    'n_communities': 0,
                    'eigenvalues': np.array([]),
                    'eigenvectors': np.array([])
                }
            
            # Convert trait-variant matrix to trait-trait matrix
            print("Converting to trait-trait network...")
            trait_trait_matrix = adj_matrix @ adj_matrix.T
            print(f"Trait-trait matrix shape: {trait_trait_matrix.shape}")
            
            # Calculate degree matrix using trait-trait matrix
            degrees = trait_trait_matrix.sum(axis=1).A1  # Convert matrix to 1D array
            m = trait_trait_matrix.sum() / 2
            
            if m == 0:
                logger.warning("No edges in network - returning default modularity values")
                return {
                    'modularity_score': 0,
                    'n_communities': 0,
                    'eigenvalues': np.array([]),
                    'eigenvectors': np.array([])
                }
            
            # Calculate expected connections
            expected = np.outer(degrees, degrees) / (2 * m)
            
            # Convert sparse matrix to dense for matrix operations
            adj_dense = trait_trait_matrix.todense()
            
            # Calculate modularity matrix
            print("Calculating modularity matrix...")
            modularity_matrix = adj_dense - expected
            
            try:
                # Perform eigendecomposition
                print("Performing eigendecomposition...")
                eigenvalues, eigenvectors = np.linalg.eigh(modularity_matrix)
                n_communities = sum(eigenvalues > 0)
                
                # Calculate overall modularity score
                Q = np.trace(modularity_matrix) / (2 * m)
                
                print(f"Analysis complete. Found {n_communities} communities with modularity score {Q}")
                
            except np.linalg.LinAlgError:
                logger.warning("Error in eigendecomposition - using fallback calculation")
                eigenvalues = np.array([])
                eigenvectors = np.array([])
                n_communities = 0
                Q = 0
            
            return {
                'modularity_score': float(Q),
                'n_communities': int(n_communities),
                'eigenvalues': eigenvalues,
                'eigenvectors': eigenvectors
            }
                
        except Exception as e:
            logger.error(f"Error analyzing network modularity: {e}")
            raise

    def _identify_population_effects(self, network: Dict, data: pd.DataFrame) -> Dict:
        """Identify population-specific effects in trait associations"""
        try:
            results = {}
            
            # Analyze each population separately
            for pop in data['population'].unique():
                pop_data = data[data['population'] == pop]
                
                # Calculate population-specific trait associations
                pop_adj = np.zeros_like(network['adjacency_matrix'])
                
                for trait_idx, trait in enumerate(network['traits']):
                    for var_idx, variant in enumerate(network['variants']):
                        mask = pop_data['variant_id'] == variant
                        if mask.any():
                            var_traits = pop_data.loc[mask, 'traits'].iloc[0]
                            if trait in var_traits:
                                pop_adj[trait_idx, var_idx] = pop_data.loc[mask, 'allele_freq'].iloc[0]
                
                # Calculate population-specific metrics
                results[pop] = {
                    'trait_connectivity': pop_adj.sum(axis=1),
                    'variant_impact': pop_adj.sum(axis=0),
                    'unique_associations': len(np.nonzero(pop_adj)[0])
                }
            
            return results
            
        except Exception as e:
            logger.error(f"Error identifying population effects: {e}")
            raise

    def _calculate_context_dependency(self, network: Dict, pop_effects: Dict) -> Dict:
        """Calculate context dependency metrics for trait associations"""
        try:
            # Initialize results
            context_scores = {}
            
            # Get baseline connectivity from network
            baseline_connectivity = network['adjacency_matrix'].sum(axis=1)
            
            # Calculate variation across populations
            pop_variations = {}
            for pop, effects in pop_effects.items():
                # Calculate relative difference from baseline
                rel_diff = (effects['trait_connectivity'] - baseline_connectivity) / baseline_connectivity
                pop_variations[pop] = {
                    'relative_difference': rel_diff,
                    'mean_deviation': np.mean(np.abs(rel_diff)),
                    'max_deviation': np.max(np.abs(rel_diff))
                }
            
            # Calculate context dependency scores
            for trait_idx, trait in enumerate(network['traits']):
                # Get trait-specific variations
                trait_vars = [pop_var['relative_difference'][trait_idx] 
                            for pop_var in pop_variations.values()]
                
                context_scores[trait] = {
                    'variance': np.var(trait_vars),
                    'range': np.ptp(trait_vars),
                    'consistency': 1 - (np.std(trait_vars) / np.mean(np.abs(trait_vars)))
                }
            
            return {
                'population_variations': pop_variations,
                'trait_context_scores': context_scores,
                'mean_context_dependency': np.mean([s['variance'] for s in context_scores.values()])
            }
            
        except Exception as e:
            logger.error(f"Error calculating context dependency: {e}")
            raise
    
    def _visualize_pleiotropy_networks(self, output_path: Path) -> None:
        """Generate visualizations for pleiotropy network analysis"""
        if not self.pleiotropy_networks:
            logger.warning("No pleiotropy analysis results to visualize")
            return
            
        try:
            # Access network data from correct location
            network_data = self.pleiotropy_networks['network_analysis']['network']
            modularity = self.pleiotropy_networks['network_analysis']['modularity']
            '''
            print("\nDebugging modularity:")
            print(type(modularity))
            print(modularity)
            
            print("\nDebugging selection_signals:")
            print("Keys available:", self.selection_signals.keys())
            print("\nContent of selection_signals:")
            for key, value in self.selection_signals.items():
                print(f"\n{key}:")
                print(value)
            '''
            # 1. Create trait network visualization
            plt.figure(figsize=(15, 10))
            
            # Convert sparse matrix to dense for visualization if needed
            adj_matrix = network_data['adjacency_matrix'].todense()
            
            # Create network graph
            G = nx.Graph()
            
            # Add nodes (traits)
            for i, trait in enumerate(network_data['traits']):
                G.add_node(trait)
                
            # Add edges based on adjacency matrix
            for i in range(len(network_data['traits'])):
                for j in range(i+1, len(network_data['traits'])):
                    if adj_matrix[i, j] > 0:
                        G.add_edge(network_data['traits'][i], 
                                network_data['traits'][j], 
                                weight=adj_matrix[i, j])
            
            # Draw network
            pos = nx.spring_layout(G)
            nx.draw(G, pos,
                    node_size=50,
                    node_color='lightblue',
                    alpha=0.6,
                    with_labels=False)
            
            # Use a more general format for the title
            
            plt.title(f'Trait Network (Modularity: {modularity["modularity_score"]:.3f}, Communities: {modularity["n_communities"]})')
            plt.savefig(output_path / 'trait_network.png', dpi=300, bbox_inches='tight')
            plt.close()
        
        
            
            # 2. Create population effects visualization
            pop_effects = self.pleiotropy_networks['population_effects']
            pop_data = pd.DataFrame({
                'Population': list(pop_effects.keys()),
                'Pleiotropic Variants': [p['pleiotropic_variants'] for p in pop_effects.values()],
                'Mean Traits per Variant': [p['mean_traits_per_variant'] for p in pop_effects.values()]
            })
            
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
            
            sns.barplot(data=pop_data,
                    x='Population',
                    y='Pleiotropic Variants',
                    ax=ax1)
            ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45)
            ax1.set_title('Number of Pleiotropic Variants by Population')
            
            sns.barplot(data=pop_data,
                    x='Population',
                    y='Mean Traits per Variant',
                    ax=ax2)
            ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45)
            ax2.set_title('Mean Traits per Variant by Population')
            
            plt.tight_layout()
            plt.savefig(output_path / 'population_pleiotropy.png', dpi=300)
            plt.close()
            
        except Exception as e:
            logger.error(f"Error visualizing pleiotropy networks: {e}")
            raise
        
    def _calculate_ihs_scores(self, data: pd.DataFrame) -> pd.DataFrame:
        """Calculate integrated haplotype scores (iHS)"""
        try:
            # Group data by population
            ihs_results = []
            
            for pop in data['population'].unique():
                pop_data = data[data['population'] == pop].copy()
                
                # Calculate derived allele frequency
                freqs = pop_data['allele_freq']
                
                # Calculate haplotype diversity
                hap_div = -freqs * np.log(freqs) - (1-freqs) * np.log(1-freqs)
                
                # Calculate iHS
                ihs = pd.DataFrame({
                    'population': pop,
                    'variant_id': pop_data['variant_id'],
                    'ihs_score': hap_div / np.std(hap_div),
                    'position': pop_data['position'],
                    'allele_freq': freqs
                })
                
                ihs_results.append(ihs)
            
            return pd.concat(ihs_results, ignore_index=True)
            
        except Exception as e:
            logger.error(f"Error calculating iHS scores: {e}")
            raise

    def _calculate_xpehh_scores(self, data: pd.DataFrame) -> pd.DataFrame:
        """Calculate cross-population extended haplotype homozygosity (XP-EHH)"""
        try:
            xpehh_results = []
            
            # Get African and non-African populations
            african_pops = data[data['major_group'] == 'African']['population'].unique()
            other_pops = data[data['major_group'] != 'African']['population'].unique()
            
            # Calculate XP-EHH for each African vs non-African population pair
            for af_pop in african_pops:
                for other_pop in other_pops:
                    # Get population data
                    pop1_data = data[data['population'] == af_pop]
                    pop2_data = data[data['population'] == other_pop]
                    
                    # Calculate EHH for each population
                    ehh1 = -np.log(pop1_data['allele_freq'])
                    ehh2 = -np.log(pop2_data['allele_freq'])
                    
                    # Calculate XP-EHH
                    xpehh = pd.DataFrame({
                        'population_1': af_pop,
                        'population_2': other_pop,
                        'variant_id': pop1_data['variant_id'],
                        'xpehh_score': np.log(ehh1/ehh2),
                        'position': pop1_data['position']
                    })
                    
                    xpehh_results.append(xpehh)
            
            return pd.concat(xpehh_results, ignore_index=True)
            
        except Exception as e:
            logger.error(f"Error calculating XP-EHH scores: {e}")
            raise

    def _analyze_environmental_correlation(self, data: pd.DataFrame) -> Dict:
        """Analyze correlation between selection signals and environmental factors"""
        try:
            # Define environmental variables (example)
            env_vars = {
                'temperature': [20, 25, 30, 35, 28],  # Example values
                'rainfall': [1000, 1200, 800, 600, 900],
                'altitude': [0, 1000, 500, 200, 300]
            }
            
            results = {}
            
            # Calculate correlations for each selection metric
            selection_metrics = {
                'ihs': data.groupby('population')['ihs_score'].mean(),
                'xpehh': data.groupby('population')['xpehh_score'].mean()
            }
            
            for metric_name, metric_values in selection_metrics.items():
                correlations = {}
                for env_var, env_values in env_vars.items():
                    corr, p_value = stats.pearsonr(metric_values, env_values)
                    correlations[env_var] = {
                        'correlation': corr,
                        'p_value': p_value
                    }
                results[metric_name] = correlations
            
            return results
            
        except Exception as e:
            logger.error(f"Error analyzing environmental correlations: {e}")
            raise

    def _identify_selective_sweeps(self, ihs_scores: pd.DataFrame, 
                                xpehh_scores: pd.DataFrame) -> Dict:
        """Identify selective sweeps using iHS and XP-EHH scores"""
        try:
            sweeps = {}
            
            # Identify sweeps using iHS
            ihs_threshold = 2.0  # Standard threshold for significance
            significant_ihs = ihs_scores[abs(ihs_scores['ihs_score']) > ihs_threshold]
            
            # Group by population
            for pop in significant_ihs['population'].unique():
                pop_sweeps = significant_ihs[significant_ihs['population'] == pop]
                
                # Identify sweep regions
                sweep_regions = []
                current_sweep = []
                
                for _, variant in pop_sweeps.iterrows():
                    if not current_sweep or \
                    variant['position'] - current_sweep[-1]['position'] < 100000:  # 100kb window
                        current_sweep.append(variant)
                    else:
                        if len(current_sweep) >= 3:  # Minimum 3 variants for a sweep
                            sweep_regions.append(current_sweep)
                        current_sweep = [variant]
                
                sweeps[pop] = {
                    'n_sweeps': len(sweep_regions),
                    'regions': [{
                        'start': min(sweep['position'] for sweep in region),
                        'end': max(sweep['position'] for sweep in region),
                        'n_variants': len(region),
                        'mean_ihs': np.mean([sweep['ihs_score'] for sweep in region])
                    } for region in sweep_regions]
                }
            
            return sweeps
            
        except Exception as e:
            logger.error(f"Error identifying selective sweeps: {e}")
            raise

    def _visualize_selection_patterns(self, output_path: Path) -> None:
        """Generate visualizations for selection pattern analysis"""
        if not self.selection_signals:
            logger.warning("No selection pattern results to visualize")
            return
            
        try:
            # Access data with correct path
            ihs_data = self.selection_signals['selection_statistics']['ihs_scores']
            xpehh_data = self.selection_signals['selection_statistics']['xpehh_scores']
            env_corr = self.selection_signals['environmental_correlations']
            summary = self.selection_signals['summary_metrics']
            
            # 1. Create iHS score distribution plot
            plt.figure(figsize=(12, 6))
            
            # Filter out NaN values for iHS scores
            valid_ihs = ihs_data.dropna(subset=['ihs_score'])
            if len(valid_ihs) > 0:
                sns.boxplot(data=valid_ihs, 
                        x='population', 
                        y='ihs_score')
                plt.xticks(rotation=45)
                plt.title('Distribution of iHS Scores by Population')
                plt.xlabel('Population')
                plt.ylabel('iHS Score')
                
                plt.tight_layout()
                plt.savefig(output_path / 'ihs_distribution.png', dpi=300)
            else:
                logger.warning("No valid iHS scores for visualization")
            plt.close()
            
            # 2. Create XP-EHH comparison plot
            plt.figure(figsize=(12, 6))
            
            # Filter out NaN values for XP-EHH scores
            valid_xpehh = xpehh_data.dropna(subset=['xpehh_score', 'variant_id'])
            if len(valid_xpehh) > 0:
                sns.boxplot(data=valid_xpehh,
                        x='population_1',
                        y='xpehh_score',
                        hue='population_2')
                plt.xticks(rotation=45)
                plt.title('XP-EHH Score Distribution')
                plt.xlabel('Population 1')
                plt.ylabel('XP-EHH Score')
                
                plt.tight_layout()
                plt.savefig(output_path / 'xpehh_comparison.png', dpi=300)
            else:
                logger.warning("No valid XP-EHH scores for visualization")
            plt.close()
            
            # 3. Create environmental correlation heatmap
            plt.figure(figsize=(10, 8))
            env_data = {}
            
            for metric, correlations in env_corr.items():
                env_data[metric] = {
                    env: corr['correlation'] 
                    for env, corr in correlations.items()
                }
            
            env_corr_df = pd.DataFrame(env_data)
            
            sns.heatmap(env_corr_df, 
                    annot=True, 
                    cmap='RdBu_r',
                    center=0,
                    fmt='.2f')
            plt.title('Environmental Correlations')
            plt.tight_layout()
            plt.savefig(output_path / 'environmental_correlations.png', dpi=300)
            plt.close()
            
            # 4. Summary metrics visualization
            plt.figure(figsize=(10, 6))
            
            metrics = {
                'Total Protective\nVariants': summary['total_protective_variants'],
                'African\nSpecific': summary['african_specific_count'],
                'Selective\nSweeps': summary['selective_sweep_count']
            }
            
            # Create bar plot
            colors = ['#2ecc71', '#3498db', '#e74c3c']  # Green, Blue, Red
            bars = plt.bar(metrics.keys(), metrics.values(), color=colors)
            
            # Customize plot
            plt.title('Selection Pattern Summary Metrics', pad=20)
            plt.ylabel('Count')
            
            # Add value labels on top of bars
            for bar in bars:
                height = bar.get_height()
                plt.text(bar.get_x() + bar.get_width()/2., height,
                        f'{int(height)}',
                        ha='center', va='bottom')
            
            plt.tight_layout()
            plt.savefig(output_path / 'selection_summary.png', dpi=300)
            plt.close()
            
            logger.info("Selection pattern visualizations completed successfully")
            
        except Exception as e:
            logger.error(f"Error visualizing selection patterns: {e}")
            raise
        
    def _random_subsample(self, data: pd.DataFrame, size: int) -> pd.DataFrame:
        """Create random subsample of specified size"""
        try:
            # Group by population to maintain population structure
            subsamples = []
            for pop in data['population'].unique():
                pop_data = data[data['population'] == pop]
                if len(pop_data) > size:
                    pop_subsample = pop_data.sample(n=size)
                else:
                    pop_subsample = pop_data  # Keep all data if population is smaller than size
                subsamples.append(pop_subsample)
                
            return pd.concat(subsamples, ignore_index=True)
            
        except Exception as e:
            logger.error(f"Error in subsampling: {e}")
            raise

    def _calculate_detection_metrics(self, data: pd.DataFrame) -> Dict:
        """Calculate variant detection metrics for a dataset"""
        try:
            metrics = {
                'total_variants': len(data['variant_id'].unique()),
                'rare_variants': len(data[data['allele_freq'] < 0.01]['variant_id'].unique()),
                'common_variants': len(data[data['allele_freq'] > 0.05]['variant_id'].unique()),
                'mean_frequency': data['allele_freq'].mean(),
                'median_frequency': data['allele_freq'].median(),
                'pathogenic_variants': len(data[data['is_pathogenic']]['variant_id'].unique())
            }
            
            # Calculate variant density
            genomic_range = data['position'].max() - data['position'].min()
            metrics['variant_density'] = metrics['total_variants'] / (genomic_range / 1e6)  # per Mb
            
            return metrics
            
        except Exception as e:
            logger.error(f"Error calculating detection metrics: {e}")
            raise

    def _calculate_detection_power(self, data: pd.DataFrame) -> Dict:
        """Calculate statistical power for variant detection"""
        try:
            # Define frequency bins
            freq_bins = [0, 0.001, 0.01, 0.05, 0.1, 1.0]
            power_results = {}
            
            # Calculate power for each frequency bin
            for i in range(len(freq_bins)-1):
                bin_data = data[(data['allele_freq'] > freq_bins[i]) & 
                            (data['allele_freq'] <= freq_bins[i+1])]
                
                # Calculate detection rate
                detected = len(bin_data[bin_data['sample_size'] >= 30])  # Minimum sample size threshold
                total = len(bin_data)
                
                power_results[f"freq_{freq_bins[i]:.3f}_{freq_bins[i+1]:.3f}"] = {
                    'detection_rate': detected/total if total > 0 else 0,
                    'n_variants': total,
                    'n_detected': detected
                }
                
            return power_results
            
        except Exception as e:
            logger.error(f"Error calculating detection power: {e}")
            raise

    def _estimate_minimum_samples(self, power_analysis: Dict) -> Dict:
        """Estimate minimum sample sizes needed for different variant frequencies"""
        try:
            min_samples = {}
            
            # Calculate for different detection thresholds
            detection_thresholds = [0.8, 0.9, 0.95]
            
            for freq_bin, stats in power_analysis.items():
                min_samples[freq_bin] = {}
                
                # Theoretical minimum sample size for each threshold
                for threshold in detection_thresholds:
                    # Using standard statistical power calculation
                    min_n = np.ceil(np.log(1-threshold) / 
                                np.log(1-float(freq_bin.split('_')[1])))
                    
                    min_samples[freq_bin][f'threshold_{threshold}'] = {
                        'min_samples': int(min_n),
                        'actual_detection_rate': stats['detection_rate']
                    }
                    
            return min_samples
            
        except Exception as e:
            logger.error(f"Error estimating minimum samples: {e}")
            raise

    def _visualize_sample_size_effects(self, output_path: Path) -> None:
        """Generate visualizations for sample size analysis"""
        if not self.sample_size_effects:
            logger.warning("No sample size analysis results to visualize")
            return
                
        try:
            print("\nAvailable keys in sample_size_effects:")
            print(self.sample_size_effects.keys())
            
            # 1. Create detection rate vs sample size plot
            plt.figure(figsize=(12, 6))
            
            if 'detection_analysis' in self.sample_size_effects:
                detection_data = self.sample_size_effects['detection_analysis']
                try:
                    # Extract subsample results and create DataFrame
                    subsample_data = []
                    for size, metrics in detection_data['population_detection'].items():
                        subsample_data.append({
                            'sample_size': metrics['sample_size'],
                            'total_variants': metrics['rare_variants_detected'],
                            'rare_variants': metrics['rare_variants_detected']
                        })
                    subsample_results = pd.DataFrame(subsample_data)
                    
                    plt.plot(subsample_results['sample_size'], 
                            subsample_results['total_variants'],
                            label='All Variants',
                            marker='o')
                    plt.plot(subsample_results['sample_size'],
                            subsample_results['rare_variants'],
                            label='Rare Variants',
                            marker='s')
                    
                    plt.title('Variant Detection vs Sample Size')
                    plt.xlabel('Sample Size')
                    plt.ylabel('Number of Variants Detected')
                    plt.legend()
                    plt.grid(True, alpha=0.3)
                    
                    plt.savefig(output_path / 'detection_rate.png', dpi=300, bbox_inches='tight')
                except Exception as e:
                    logger.warning(f"Error creating detection rate plot: {e}")
            else:
                logger.warning("No detection analysis data available")
            plt.close()
            
            # 2. Create power analysis heatmap if data available
            if 'power_analysis' in self.sample_size_effects:
                plt.figure(figsize=(10, 8))
                power_data = self.sample_size_effects['power_analysis']
                
                try:
                    # Extract power curves data
                    power_curves = power_data['power_curves']
                    sample_sizes = sorted(set(curve[0]['sample_size'] for curve in power_curves.values()))
                    freq_bins = list(power_curves.keys())
                    
                    # Create power matrix
                    power_matrix = np.zeros((len(sample_sizes), len(freq_bins)))
                    for i, freq_bin in enumerate(freq_bins):
                        curve = power_curves[freq_bin]
                        for j, size in enumerate(sample_sizes):
                            size_data = next((d for d in curve if d['sample_size'] == size), None)
                            power_matrix[j, i] = size_data['power'] if size_data else np.nan
                    
                    power_df = pd.DataFrame(power_matrix, 
                                        index=sample_sizes, 
                                        columns=[bin.replace('_', '-') for bin in freq_bins])
                    
                    sns.heatmap(power_df,
                            annot=True,
                            fmt='.2f',
                            cmap='YlOrRd',
                            cbar_kws={'label': 'Detection Power'})
                    
                    plt.title('Detection Power by Sample Size and Frequency')
                    plt.xlabel('Allele Frequency Range')
                    plt.ylabel('Sample Size')
                    plt.tight_layout()
                    
                    plt.savefig(output_path / 'power_analysis.png', dpi=300, bbox_inches='tight')
                except Exception as e:
                    logger.warning(f"Error creating power analysis heatmap: {e}")
            plt.close()
            
            # 3. Create minimum sample size recommendations plot
            if 'minimum_requirements' in self.sample_size_effects:
                plt.figure(figsize=(12, 6))
                
                try:
                    min_samples = self.sample_size_effects['minimum_requirements']
                    sample_data = []
                    
                    for freq_bin, thresholds in min_samples.items():
                        for threshold, stats in thresholds.items():
                            if isinstance(stats, dict) and 'min_samples' in stats:
                                sample_data.append({
                                    'frequency_bin': freq_bin.replace('_', '-'),
                                    'threshold': float(threshold.split('_')[1]) if '_' in threshold else 0,
                                    'min_samples': stats['min_samples']
                                })
                    
                    if sample_data:
                        sample_df = pd.DataFrame(sample_data)
                        
                        # Create grouped bar plot
                        g = sns.barplot(data=sample_df,
                                    x='frequency_bin',
                                    y='min_samples',
                                    hue='threshold',
                                    palette='viridis')
                        
                        plt.title('Minimum Sample Size Requirements by Frequency')
                        plt.xlabel('Allele Frequency Range')
                        plt.ylabel('Minimum Sample Size')
                        plt.xticks(rotation=45)
                        
                        # Add value labels on bars
                        for container in g.containers:
                            g.bar_label(container, fmt='%d')
                        
                        plt.legend(title='Power Threshold', bbox_to_anchor=(1.05, 1), loc='upper left')
                        plt.tight_layout()
                        
                        plt.savefig(output_path / 'minimum_samples.png', dpi=300, bbox_inches='tight')
                    else:
                        logger.warning("No valid minimum sample size data")
                except Exception as e:
                    logger.warning(f"Error creating minimum sample size plot: {e}")
            plt.close()
            
            # 4. Create population-specific detection rates plot
            if 'rare_variant_assessment' in self.sample_size_effects:
                plt.figure(figsize=(12, 6))
                
                try:
                    pop_detection = self.sample_size_effects['rare_variant_assessment']['population_detection']
                    pop_data = pd.DataFrame([
                        {
                            'Population': pop,
                            'Detection Rate': stats['detection_rate'],
                            'Sample Size': stats['sample_size']
                        }
                        for pop, stats in pop_detection.items()
                    ])
                    
                    # Create scatter plot
                    sns.scatterplot(data=pop_data,
                                x='Sample Size',
                                y='Detection Rate',
                                s=100)
                    
                    # Add population labels
                    for _, row in pop_data.iterrows():
                        plt.annotate(row['Population'],
                                (row['Sample Size'], row['Detection Rate']),
                                xytext=(5, 5),
                                textcoords='offset points')
                    
                    plt.title('Population-Specific Rare Variant Detection Rates')
                    plt.xlabel('Sample Size')
                    plt.ylabel('Detection Rate')
                    plt.grid(True, alpha=0.3)
                    plt.tight_layout()
                    
                    plt.savefig(output_path / 'population_detection.png', dpi=300, bbox_inches='tight')
                except Exception as e:
                    logger.warning(f"Error creating population detection plot: {e}")
            plt.close()
            
            logger.info("Sample size effect visualizations completed successfully")
                
        except Exception as e:
            logger.error(f"Error visualizing sample size effects: {e}")
            raise
        

        
    def save_analysis_results(extended_analysis, filename='analysis_results.pkl'):
        """Save all analysis results to a file"""
        results_to_save = {
            'geo_clusters': extended_analysis.geo_clusters,
            'pleiotropy_networks': extended_analysis.pleiotropy_networks,
            'selection_signals': extended_analysis.selection_signals,
            'sample_size_effects': extended_analysis.sample_size_effects
        }
        with open(filename, 'wb') as f:
            pickle.dump(results_to_save, f)
        print(f"Analysis results saved to {filename}")

    def load_analysis_results(self, filename: str) -> Dict:
        """
        Load analysis results from file and setup both object attributes and return results dict
        
        Args:
            filename: Path to the analysis results PKL file
            
        Returns:
            Dict containing all analysis results
        """
        try:
            if not os.path.exists(filename):
                raise FileNotFoundError(f"Analysis results file {filename} not found. Run analyze mode first.")
            
            with open(filename, 'rb') as f:
                results = pickle.load(f)
                
            # Add detailed debugging information
            logger.info("Type of loaded results: %s", type(results))
            logger.info("All available keys in results: %s", results.keys())
            
            for key, value in results.items():
                logger.info(f"\nKey: {key}")
                logger.info(f"Type: {type(value)}")
                if isinstance(value, dict):
                    logger.info(f"Nested keys: {value.keys()}")
                if hasattr(value, 'shape'):
                    logger.info(f"Shape: {value.shape}")
                if hasattr(value, 'columns'):
                    logger.info(f"Columns: {value.columns.tolist()}")
                # If it's a DataFrame, show first few rows
                if hasattr(value, 'head'):
                    logger.info(f"First few rows:\n{value.head()}")
            
            # Load results into object attributes
            self.geo_clusters = results.get('geo_clusters', {})
            self.pleiotropy_networks = results.get('pleiotropy_networks', {})
            self.selection_signals = results.get('selection_signals', {})
            self.sample_size_effects = results.get('sample_size_effects', {})
            
            # Create results dictionary
            report_results = {
                'geo_clusters': self.geo_clusters,
                'pleiotropy_networks': self.pleiotropy_networks,
                'selection_signals': self.selection_signals,
                'sample_size_effects': self.sample_size_effects
            }
            
            logger.info("Loaded results keys: %s", list(report_results.keys()))
            logger.info("Analysis results loaded from %s", filename)
            
            return report_results
                
        except Exception as e:
            logger.error(f"Error loading analysis results: {e}")
            logger.error("Full error details:", exc_info=True)
            raise
    def _write_geographic_clusters_section(self, f: TextIO, geo_results: Dict) -> None:
        """Write geographic clustering analysis following schema structure"""
        try:
            # Research Question
            f.write("## Geographic Clustering Analysis\n\n")
            f.write("### Research Question\n")
            f.write("How do risk alleles cluster geographically across African populations, and what ")
            f.write("is their relationship with documented disease prevalence patterns?\n\n")

            # Key Metrics/Statistics
            stats = geo_results.get('spatial_statistics', {})
            f.write("### Key Metrics\n")
            
            # Moran's I statistics
            moran_i = stats.get('global_moran_i')
            if isinstance(moran_i, (float, np.float64)):
                f.write(f"- Moran's I: {float(moran_i):.4f}\n")
            else:
                f.write("- Moran's I: Not calculated\n")

            # P-value
            pvalue = stats.get('significance', {}).get('p_value')
            if isinstance(pvalue, (float, np.float64)):
                f.write(f"- Statistical significance: p = {float(pvalue):.4e}\n")
            else:
                f.write("- Statistical significance: Not calculated\n")

            # Cluster Analysis
            cluster_analysis = geo_results.get('cluster_analysis', {})
            n_clusters = cluster_analysis.get('number_of_clusters', 0)
            f.write(f"- Number of clusters: {n_clusters}\n")
            
            mean_size = cluster_analysis.get('mean_cluster_size')
            if isinstance(mean_size, (float, np.float64)):
                f.write(f"- Mean cluster size: {float(mean_size):.2f}\n\n")
            else:
                f.write("- Mean cluster size: Not calculated\n\n")

            # Analysis Results
            f.write("### Analysis Results\n")

            # Spatial Clusters
            f.write("#### Spatial Clustering Patterns\n")
            clusters_list = cluster_analysis.get('clusters', [])
            for cluster in clusters_list:
                if isinstance(cluster, dict):
                    f.write(f"Cluster in {cluster.get('region', 'Unknown')}:\n")
                    f.write(f"- Type: {cluster.get('type', 'Unknown')}\n")
                    mean_freq = cluster.get('mean_frequency')
                    if isinstance(mean_freq, (float, np.float64)):
                        f.write(f"- Mean frequency: {float(mean_freq):.4f}\n")
                    else:
                        f.write("- Mean frequency: Not calculated\n")
                    f.write(f"- Variants: {cluster.get('n_variants', 0)}\n\n")

            # Environmental Correlations
            f.write("#### Environmental Correlations\n")
            env_corr = geo_results.get('environmental_correlations', {})
            for env_var, corr in env_corr.items():
                if isinstance(corr, dict):
                    f.write(f"{env_var}:\n")
                    correlation = corr.get('correlation')
                    if isinstance(correlation, (float, np.float64)):
                        f.write(f"- Correlation: {float(correlation):.4f}\n")
                    pval = corr.get('p_value')
                    if isinstance(pval, (float, np.float64)):
                        f.write(f"- Significance: p = {float(pval):.4e}\n")
                    f.write("\n")
            
            # Disease Prevalence Relationships
            f.write("#### Disease Prevalence Relationships\n")
            correlations = geo_results.get('disease_correlations', {})
            for disease, stats in correlations.items():
                f.write(f"\n{disease}:\n")
                if isinstance(stats, dict):
                    corr = stats.get('correlation')
                    if isinstance(corr, (float, np.float64)):
                        f.write(f"- Correlation: {float(corr):.4f}\n")
                    pval = stats.get('p_value')
                    if isinstance(pval, (float, np.float64)):
                        f.write(f"- Significance: p = {float(pval):.4e}\n")
            f.write("\n")

            # Population-Specific Findings
            f.write("### Population-Specific Findings\n")
            pop_patterns = geo_results.get('population_patterns', {})
            
            # Regional Distribution
            f.write("#### Regional Distribution Patterns\n")
            for region, patterns in pop_patterns.items():
                if isinstance(patterns, dict):
                    f.write(f"\n{region}:\n")
                    mean_freq = patterns.get('mean_frequency')
                    if isinstance(mean_freq, (float, np.float64)):
                        f.write(f"- Mean frequency: {float(mean_freq):.4f}\n")
                    f.write(f"- Variant count: {patterns.get('variant_count', 0)}\n")
                    f.write(f"- High-risk variants: {patterns.get('high_risk_variants', 0)}\n")

            # Population-Specific Variants
            f.write("\n#### Population-Specific Variants\n")
            for pop, patterns in pop_patterns.items():
                if isinstance(patterns, dict) and patterns.get('unique_variants', 0) > 0:
                    f.write(f"{pop}: {patterns['unique_variants']} unique variants\n")

        except Exception as e:
            logger.error(f"Error writing geographic clusters section: {e}")
            logger.error(f"Error details: {traceback.format_exc()}")
            f.write("\nError occurred while writing geographic clusters section.\n")
    
    def _write_pleiotropy_section(self, f: TextIO, pleio_results: Dict) -> None:
        """Write pleiotropic effects analysis following schema structure"""
        try:
            # Research Question
            f.write("## Pleiotropic Effects Analysis\n\n")
            f.write("### Research Question\n")
            f.write("What does the distribution of pleiotropic variants reveal about the ")
            f.write("context-dependency of genetic effects in different populations?\n\n")

            # Key Metrics/Statistics
            f.write("### Key Metrics\n")
            network = pleio_results.get('network_analysis', {})
            
            # Network metrics
            modularity = network.get('modularity')
            if isinstance(modularity, (float, np.float64)):
                f.write(f"- Network modularity: {float(modularity):.4f}\n")
            else:
                f.write("- Network modularity: Not calculated\n")
                
            size = len(network.get('traits', [])) if 'traits' in network else 0
            f.write(f"- Network size: {size} nodes\n")
            
            density = network.get('density', 0)
            if isinstance(density, (float, np.float64)):
                f.write(f"- Network density: {float(density):.4f}\n")

            # Trait Connection Statistics
            trait_stats = network.get('trait_connection_stats', {})
            f.write("\n#### Trait Connection Statistics\n")
            if trait_stats:
                f.write(f"- Mean connections per trait: {trait_stats.get('mean_connections', 0):.2f}\n")
                f.write(f"- Maximum connections: {trait_stats.get('max_connections', 0)}\n")
                f.write(f"- Total connections: {trait_stats.get('total_connections', 0)}\n")

            # Statistical Tests
            f.write("\n#### Statistical Significance Tests\n")
            tests = pleio_results.get('statistical_tests', {})
            if 'context_dependency' in tests:
                test = tests['context_dependency']
                f.write(f"- Test Method: {test.get('method', 'Not specified')}\n")
                stat = test.get('statistic')
                if isinstance(stat, (float, np.float64)):
                    f.write(f"- Test Statistic: {float(stat):.4f}\n")
                pval = test.get('p_value')
                if isinstance(pval, (float, np.float64)):
                    f.write(f"- P-value: {float(pval):.4e}\n")
            f.write("\n")

            # Analysis Results
            f.write("### Analysis Results\n")
            
            # Trait Network Characteristics
            f.write("#### Trait Network Characteristics\n")
            if 'traits' in network:
                f.write(f"- Total traits: {len(network['traits'])}\n")
                f.write(f"- Total variants: {len(network.get('variants', []))}\n")

            # Context Dependency Scores
            f.write("\n#### Context Dependency Scores\n")
            context_scores = pleio_results.get('context_scores', {})
            for pop, scores in context_scores.items():
                if isinstance(scores, dict):
                    f.write(f"{pop}:\n")
                    mean_dev = scores.get('mean_deviation')
                    if isinstance(mean_dev, (float, np.float64)):
                        f.write(f"- Mean deviation: {float(mean_dev):.4f}\n")
                    cons_score = scores.get('consistency_score')
                    if isinstance(cons_score, (float, np.float64)):
                        f.write(f"- Consistency score: {float(cons_score):.4f}\n")

            # Network Community Structure
            f.write("\n#### Network Community Structure\n")
            communities = network.get('communities', [])
            for i, comm in enumerate(communities, 1):
                if isinstance(comm, dict):
                    f.write(f"Community {i}:\n")
                    f.write(f"- Size: {comm.get('size', 0)} nodes\n")
                    if 'density' in comm:
                        f.write(f"- Density: {float(comm['density']):.4f}\n")

            # Population-Specific Findings
            f.write("\n### Population-Specific Findings\n")
            
            # Population Effects
            f.write("#### Population-Specific Trait Associations\n")
            pop_effects = pleio_results.get('population_effects', {})
            for pop, effects in pop_effects.items():
                if isinstance(effects, dict):
                    f.write(f"\n{pop}:\n")
                    f.write(f"- Unique traits: {effects.get('unique_traits', 0)}\n")
                    f.write(f"- Maximum traits per variant: {effects.get('max_traits_per_variant', 0)}\n")
                    f.write(f"- Pleiotropic variants: {effects.get('pleiotropic_variants', 0)}\n")

            # Pleiotropic Patterns
            f.write("\n#### Unique Pleiotropic Patterns\n")
            sig_diffs = pleio_results.get('statistical_tests', {}).get('significant_differences', [])
            for diff in sig_diffs:
                if isinstance(diff, dict):
                    trait_pair = diff.get('trait_pair', ('Unknown', 'Unknown'))
                    f.write(f"Trait pair: {' - '.join(trait_pair)}\n")
                    af_mean = diff.get('african_mean')
                    if isinstance(af_mean, (float, np.float64)):
                        f.write(f"- African mean: {float(af_mean):.4f}\n")
                    other_mean = diff.get('other_mean')
                    if isinstance(other_mean, (float, np.float64)):
                        f.write(f"- Other mean: {float(other_mean):.4f}\n")
                    pval = diff.get('p_value')
                    if isinstance(pval, (float, np.float64)):
                        f.write(f"- Significance: p = {float(pval):.4e}\n")
                    f.write("\n")

        except Exception as e:
            logger.error(f"Error writing pleiotropy section: {e}")
            logger.error(f"Error details: {traceback.format_exc()}")
            f.write("\nError occurred while writing pleiotropy section.\n")
    def _write_selection_patterns_section(self, f: TextIO, sel_results: Dict) -> None:
        """Write selection patterns analysis following schema structure"""
        try:
            # Research Question
            f.write("## Selection Patterns Analysis\n\n")
            f.write("### Research Question\n")
            f.write("How are protective variants distributed across populations, and what insights ")
            f.write("do these patterns provide about historical selective pressures?\n\n")

            # Key Metrics/Statistics
            f.write("### Key Metrics\n")
            stats = sel_results.get('selection_statistics', {})
            
            # iHS Statistics
            ihs = stats.get('ihs_scores', {})
            f.write("#### iHS Scores\n")
            f.write(f"- Mean iHS: {ihs.get('mean_score', 0):.4f}\n")
            f.write(f"- Significant signals: {ihs.get('significant_count', 0)}\n")

            # XP-EHH Statistics
            xpehh = stats.get('xpehh_scores', {})
            f.write("\n#### XP-EHH Scores\n")
            f.write(f"- Mean XP-EHH: {xpehh.get('mean_score', 0):.4f}\n")
            f.write(f"- Cross-population signals: {xpehh.get('significant_count', 0)}\n\n")

            # Analysis Results
            f.write("### Analysis Results\n")

            # Selection Signals
            f.write("#### Selection Signal Distribution\n")
            sweeps = sel_results.get('selective_sweeps', {})
            for pop, sweep in sweeps.items():
                f.write(f"\n{pop}:\n")
                f.write(f"- High frequency sweeps: {sweep.get('high_freq_protective', 0)}\n")
                f.write(f"- Mean frequency: {sweep.get('mean_protective_freq', 0):.4f}\n")

            # Environmental Correlations
            f.write("\n#### Environmental Correlations\n")
            env_corr = sel_results.get('environmental_correlations', {})
            for metric, correlations in env_corr.items():
                f.write(f"\n{metric} correlations:\n")
                for env_var, stats in correlations.items():
                    f.write(f"- {env_var}:\n")
                    f.write(f"  * Correlation: {stats.get('correlation', 0):.4f}\n")
                    f.write(f"  * Significance: p = {stats.get('p_value', 1):.4e}\n")

            # Population Comparisons
            f.write("\n#### Population Comparisons\n")
            comparisons = sel_results.get('population_comparisons', {})
            for test_name, test in comparisons.items():
                f.write(f"\n{test_name}:\n")
                f.write(f"- Method: {test.get('method', 'Not specified')}\n")
                f.write(f"- Statistic: {test.get('statistic', 0):.4f}\n")
                f.write(f"- Significance: p = {test.get('p_value', 1):.4e}\n")

            # Population-Specific Findings
            f.write("\n### Population-Specific Findings\n")

            # Selection Patterns
            f.write("#### Population-Specific Selection Signals\n")
            patterns = sel_results.get('selection_analysis', {}).get('patterns', {})
            for pop, data in patterns.items():
                f.write(f"\n{pop}:\n")
                if isinstance(data, dict):
                    f.write(f"- High frequency protective: {data.get('high_freq_protective', 0)}\n")
                    f.write(f"- Mean protective frequency: {data.get('mean_protective_freq', 0):.4f}\n")
                    if 'protective_consequence_types' in data:
                        f.write("- Consequence types:\n")
                        for cons, count in data['protective_consequence_types'].items():
                            f.write(f"  * {cons}: {count}\n")

            # Regional Patterns
            f.write("\n#### Regional Adaptation Patterns\n")
            for pop, data in patterns.items():
                if isinstance(data, dict) and 'adaptive_benefits' in data:
                    f.write(f"\n{pop} adaptations:\n")
                    for benefit in data['adaptive_benefits'].values():
                        f.write(f"- {benefit}\n")

        except Exception as e:
            logger.error(f"Error writing selection patterns section: {e}")
            f.write("\nError occurred while writing selection patterns section.\n")
    def _write_sample_size_effects_section(self, f: TextIO, size_results: Dict) -> None:
        """Write sample size effects analysis following schema structure"""
        try:
            # Research Question
            f.write("## Sample Size Effects Analysis\n\n")
            f.write("### Research Question\n")
            f.write("How does sample size variation across populations affect our ability to ")
            f.write("detect and characterize rare clinically significant variants?\n\n")

            # Key Metrics/Statistics
            power = size_results.get('power_analysis', {}) or {}
            if power:
                for freq_bin, stats in power.items():
                    if isinstance(stats, dict):
                        f.write(f"\n{freq_bin} frequency bin:\n")
                        detection = stats.get('detection_rate')
                        if detection is not None:
                            f.write(f"- Detection rate: {detection:.4f if isinstance(detection, (float, np.float64)) else 'N/A'}\n")
                        variants = stats.get('n_variants')
                        if variants is not None:
                            f.write(f"- Variants: {int(variants) if isinstance(variants, (int, np.integer)) else 'N/A'}\n")

            # Statistical Tests
            f.write("\n#### Statistical Significance\n")
            stats = size_results.get('population_comparisons', {}) or {}
            if stats.get('detection_rate_test'):
                test = stats['detection_rate_test']
                f.write(f"- Method: {test.get('method', 'Not specified')}\n")
                statistic = test.get('statistic')
                p_value = test.get('p_value')
                f.write(f"- Statistic: {statistic:.4f if statistic is not None else 'N/A'}\n")
                f.write(f"- Significance: p = {p_value:.4e if p_value is not None else 'N/A'}\n\n")

            # Analysis Results
            f.write("### Analysis Results\n")

            # Power Analysis Results
            f.write("#### Power Analysis Findings\n")
            power_curves = power.get('power_curves', {}) or {}
            for freq_range, curve in power_curves.items():
                if curve:  # Check if curve exists and is not empty
                    f.write(f"\n{freq_range}:\n")
                    try:
                        max_power = max(p.get('power', 0) for p in curve if isinstance(p, dict))
                        f.write(f"- Maximum power: {max_power:.4f if max_power else 'N/A'}\n")
                        
                        # Safely find sample size for 80% power
                        power_threshold = next((p.get('sample_size') for p in curve 
                                            if isinstance(p, dict) and p.get('power', 0) >= 0.8), 
                                            'Not achieved')
                        f.write(f"- Sample size for 80% power: {power_threshold}\n")
                    except Exception as e:
                        f.write(f"- Error calculating power metrics: {str(e)}\n")

            # Subsampling Results
            f.write("\n#### Subsampling Analysis\n")
            subsample = size_results.get('detection_analysis', {}).get('subsampling_results', {}) or {}
            for size, results in subsample.items():
                if isinstance(results, dict):
                    f.write(f"\nSample size {size}:\n")
                    total_variants = results.get('total_variants')
                    rare_variants = results.get('rare_variants')
                    detection_rate = results.get('detection_rate')
                    
                    f.write(f"- Total variants: {int(total_variants) if total_variants is not None else 'N/A'}\n")
                    f.write(f"- Rare variants: {int(rare_variants) if rare_variants is not None else 'N/A'}\n")
                    f.write(f"- Detection rate: {detection_rate:.4f if detection_rate is not None else 'N/A'}\n")

            # Detection Patterns
            f.write("\n#### Detection Rate Patterns\n")
            rare = size_results.get('rare_variant_assessment', {}) or {}
            total_rare = rare.get('total_rare_variants')
            mean_detection = rare.get('mean_detection_rate')
            f.write(f"- Total rare variants: {int(total_rare) if total_rare is not None else 'N/A'}\n")
            f.write(f"- Mean detection rate: {mean_detection:.4f if mean_detection is not None else 'N/A'}\n")

            # Population-Specific Findings
            f.write("\n### Population-Specific Findings\n")

            # Population Requirements
            f.write("#### Population-Specific Requirements\n")
            pop_detect = rare.get('population_detection', {}) or {}
            for pop, stats in pop_detect.items():
                if isinstance(stats, dict):
                    f.write(f"\n{pop}:\n")
                    sample_size = stats.get('sample_size')
                    variants_detected = stats.get('rare_variants_detected')
                    detect_rate = stats.get('detection_rate')
                    mean_freq = stats.get('mean_frequency')
                    
                    f.write(f"- Sample size: {sample_size:.1f if sample_size is not None else 'N/A'}\n")
                    f.write(f"- Rare variants detected: {int(variants_detected) if variants_detected is not None else 'N/A'}\n")
                    f.write(f"- Detection rate: {detect_rate:.4f if detect_rate is not None else 'N/A'}\n")
                    f.write(f"- Mean frequency: {mean_freq:.4e if mean_freq is not None else 'N/A'}\n")

            # Regional Patterns
            f.write("\n#### Regional Detection Patterns\n")
            min_reqs = size_results.get('minimum_requirements', {}) or {}
            for region, reqs in min_reqs.items():
                if isinstance(reqs, dict):
                    f.write(f"\n{region}:\n")
                    for threshold, data in reqs.items():
                        if isinstance(data, dict):
                            f.write(f"- {threshold}:\n")
                            min_samples = data.get('min_samples')
                            achieved_power = data.get('achieved_power')
                            f.write(f"  * Minimum samples: {min_samples if min_samples is not None else 'Not calculated'}\n")
                            f.write(f"  * Achieved power: {achieved_power:.4f if achieved_power is not None else 'N/A'}\n")

        except Exception as e:
            logger.error(f"Error writing sample size effects section: {e}")
            f.write("\nError occurred while writing sample size effects section.\n")
    def generate_analysis_report(self, results: Dict) -> None:
        """Generate comprehensive analysis report following schema structure"""
        try:
            report_path = Path(self.output_dir) / 'extended_analysis_report.md'
            # Debug the results
            logger.info("Available results keys: %s", list(results.keys()))
            for key, value in results.items():
                logger.info("Content summary for %s: %s", key, str(value)[:100])
                
            report_path = Path(self.output_dir) / 'extended_analysis_report.md'
            with open(report_path, 'w', encoding='utf-8') as f:
                # Title and Date
                f.write("# African Genetic Diversity Extended Analysis Report\n\n")
                f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d')}\n\n")

                # Write each analysis section in order
                self._write_geographic_clusters_section(f, results.get('geo_clusters', {}))
                f.write("\n")  # Section spacing
                
                self._write_pleiotropy_section(f, results.get('pleiotropy_networks', {}))
                f.write("\n")
                
                self._write_selection_patterns_section(f, results.get('selection_signals', {}))
                f.write("\n")
                
                self._write_sample_size_effects_section(f, results.get('sample_size_effects', {}))
                
                # Flush to ensure writing is complete
                f.flush()
                
            logger.info(f"Analysis report generated successfully: {report_path}")
            
        except Exception as e:
            logger.error(f"Error generating analysis report: {e}")
            logger.error(f"Error details: {traceback.format_exc()}")
            raise
if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(description='Extended African Genetic Diversity Analysis')
    parser.add_argument('--email', type=str, default="your@email.com", 
                        help='Your email address for API access')
    parser.add_argument('--output', type=str, default="extended_results",
                        help='Output directory for results')
    parser.add_argument('--mode', type=str, 
                        choices=['analyze', 'visualize', 'report', 'all'],
                        default='all',
                        help='Operation mode: analyze data, create visualizations, generate report, or all')
    parser.add_argument('--input', type=str, required=True,
                        help='Path to processed variant data CSV from base analysis')
    parser.add_argument('--geo_data', type=str, required=True,
                        help='Path to geographic data file')

    args = parser.parse_args()
    
    try:
        # Initialize base analysis
        logger.info("Initializing base analysis...")
        base_analysis = AfricanGeneticDiversity(
            email=args.email,
            output_dir=args.output
        )
        
        # Load existing data
        logger.info("Loading existing variant data...")
        if args.input.endswith('.csv'):
            base_analysis.variant_data = pd.read_csv(args.input, low_memory=False)
        elif args.input.endswith('.pkl'):
            logger.info(f"Loading base data from {args.input}")
            base_analysis.variant_data = pd.read_pickle(args.input)
        else:
            raise ValueError("Input file must be either .csv or .pkl")
        
        # Initialize extended analysis
        logger.info("Initializing extended analysis...")
        extended_analysis = ExtendedAfricanDiversity(
            base_analysis=base_analysis,
            geo_data_path=args.geo_data,
            output_dir=args.output
        )
        
        results = {}
        
        if args.mode in ['analyze', 'all']:
            # Run analyses and collect results
            logger.info("Running geographic clustering analysis...")
            results['geo_clusters'] = extended_analysis.analyze_geographic_clusters()
            
            logger.info("Running pleiotropy context analysis...")
            results['pleiotropy_networks'] = extended_analysis.analyze_pleiotropy_context()
            
            logger.info("Running selection patterns analysis...")
            results['selection_signals'] = extended_analysis.analyze_selection_patterns()
            
            logger.info("Running sample size effects analysis...")
            results['sample_size_effects'] = extended_analysis.analyze_sample_size_effects()
            
            # Save results after analysis
            logger.info("Saving analysis results...")
            extended_analysis.save_analysis_results(os.path.join(args.output, 'analysis_results.pkl'))
        
        # For report or visualization modes, load existing results if needed
        elif args.mode in ['report', 'visualize']:
            results_path = os.path.join(args.output, 'analysis_results.pkl')
            logger.info(f"Loading existing results from {results_path}")
            try:
                results = extended_analysis.load_analysis_results(results_path)
                if not results:
                    raise ValueError("No analysis results loaded")
            except Exception as e:
                logger.error(f"Failed to load analysis results: {e}")
                raise
        
        # Generate visualizations if requested
        if args.mode in ['visualize', 'all']:
            logger.info("Generating visualizations...")
            extended_analysis._visualize_geographic_clusters(Path(args.output))
            extended_analysis._visualize_pleiotropy_networks(Path(args.output))
            extended_analysis._visualize_selection_patterns(Path(args.output))
            extended_analysis._visualize_sample_size_effects(Path(args.output))
        
        # Generate report if requested
        if args.mode in ['report', 'all']:
            if not results:
                logger.error("No results available for report generation")
                raise ValueError("Analysis results required for report generation")
            
            logger.info("Generating analysis report...")
            extended_analysis.generate_analysis_report(results)
            logger.info("Report generation complete")
        
        logger.info(f"Operations complete. Results saved in {args.output}")
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        logger.error(f"Error details: {traceback.format_exc()}")
        raise
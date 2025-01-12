# African Genetic Diversity Analysis Pipeline

## Project Overview

This bioinformatics pipeline analyzes African genetic diversity patterns through two major analytical approaches:

Base Analysis (publication1_1.py):
├── Focuses on fundamental genetic diversity patterns across populations:
│   ├── Clinical variant impact assessment
│   ├── Functional consequence analysis
│   ├── Basic geographical distribution
│   ├── Initial pleiotropic effects
│   ├── Protective variant distribution
│   └── Sample size implications
│
└── Key Findings:
    ├── 75,216 variants analyzed across 3 population groups
    ├── 2,039 African-specific variants identified
    ├── Significant differences in variant distributions (p = 1.16e-81)
    └── 3,299 high severity variants (2,242 in African populations)

Extended Analysis (african_diversity_extended.py):
├── Builds upon base analysis with:
│   ├── Advanced geographical clustering analysis (Moran's I: 0.8454)
│   ├── Environmental correlation studies
│   ├── Detailed network-based pleiotropic analysis
│   ├── Comprehensive selection pattern analysis
│   └── Enhanced sample size effects with power analysis
│
└── Key Findings:
    ├── 9 distinct geographical clusters identified
    ├── 2,315,275 rare variants analyzed
    ├── Comprehensive power analysis across frequency bins
    └── Environmental correlation patterns

## Publications

This pipeline supports three publications:

Publication 1 (From Base Analysis):
├── Focus: Clinical impact and functional consequences
├── Key sections: Clinical Impact Analysis, Functional Consequences
└── Major finding: Population-specific variant distributions

Publication 2 (From Base Analysis):
├── Focus: Pleiotropic effects and protective variants
├── Key sections: Pleiotropic Effects, Protective Variants
└── Major finding: Population-specific trait associations

Publication 3 (From Extended Analysis):
├── Focus: Geographical patterns and selection
├── Key sections: Geographic Clustering, Selection Patterns
└── Major finding: Environmental adaptation signatures

## Repository Structure

african_diversity_project/
├── Scripts/
│   ├── Base Analysis/
│   │   ├── publication1_1.py        # Core analysis script
│   │   └── variant_cache.json       # Cached variant data
│   └── Extended Analysis/
│       ├── african_diversity_extended.py
│       └── geodata.py               # Geographic utilities
│
├── Data/
│   ├── Raw/
│   │   ├── raw_variant_data.csv     # Primary variant data
│   │   └── geo_data.csv            # Geographic reference data
│   └── Processed/
│       ├── base_analysis.pkl
│       └── extended_analysis.pkl
│
├── Reports/
│   ├── Base/
│   │   └── base_analysis_report.md
│   └── Extended/
│       └── extended_analysis_report.md
│
└── Manuscripts/
    ├── Publication1/
    ├── Publication2/
    └── Publication3/

## Key Features

Base Analysis Capabilities:
├── Population variant distribution analysis
├── Clinical significance assessment
├── Functional impact evaluation
├── Basic geographical distribution
├── Preliminary pleiotropic analysis
└── Sample size effect analysis

Extended Analysis Capabilities:
├── Spatial autocorrelation (Moran's I)
├── Environmental factor correlation
├── Network-based pleiotropy analysis
├── Selection pattern detection
└── Power analysis by frequency bin

## Usage

Base Analysis:
└── python publication1_1.py \
    --input raw_variant_data.csv \
    --output base_results

Extended Analysis:
└── python african_diversity_extended.py \
    --mode analyze \
    --input raw_variant_data.csv \
    --geo_data geo_data.csv \
    --output extended_results

## Requirements

Required Packages:
├── pandas>=1.5.0
├── numpy>=1.24.0
├── scipy>=1.9.0
├── matplotlib>=3.5.0
├── networkx>=2.8.0    # For network analysis
└── statsmodels>=0.13.0

## Data Sources

Data Sources:
├── Variant data: [source details]
├── Geographic data: [source details]
└── Population reference data: [source details]

## Result Interpretation

Analysis Reports Include:
├── Statistical summaries
├── Population-specific findings
├── Visualization outputs
└── Clinical implications

## Contributing

Contributing Guidelines:
├── Fork the repository
├── Create feature branch
├── Commit changes
├── Push to branch
└── Create Pull Request

## License

License Information:
└── [Specify license details]

## Contact

Contact Information:
└── [Your contact details]

## Acknowledgments

Acknowledgments:
├── [Funding sources]
├── [Collaborators]
└── [Supporting institutions]
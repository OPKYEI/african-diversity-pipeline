# Import both classes
from publication1_1 import AfricanGeneticDiversity
from african_diversity_extended import ExtendedAfricanDiversity

# Initialize base analysis
base_analysis = AfricanGeneticDiversity(email="your@email.com")
base_analysis.fetch_variant_data()
base_analysis.preprocess_data()

# Initialize extended analysis
extended_analysis = ExtendedAfricanDiversity(
    base_analysis=base_analysis,
    geo_data_path="path/to/geographic_data.geojson"
)

# Run new analyses
geo_results = extended_analysis.analyze_geographic_clusters()
pleio_results = extended_analysis.analyze_pleiotropy_context()
selection_results = extended_analysis.analyze_selection_patterns()
sample_results = extended_analysis.analyze_sample_size_effects()

# Generate comprehensive report
extended_analysis.generate_report()
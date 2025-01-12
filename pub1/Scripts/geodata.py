
import pandas as pd
import numpy as np

# Create a dictionary with the coordinates for each population
population_coords = {
    'ACB': {'latitude': 13.1939, 'longitude': -59.5432},  # Barbados
    'ASW': {'latitude': 39.8283, 'longitude': -98.5795},  # USA (geographic center)
    'CEU': {'latitude': 50.0755, 'longitude': 14.4378},   # Central Europe (using Prague as center)
    'CHB': {'latitude': 39.9042, 'longitude': 116.4074},  # Beijing
    'ESN': {'latitude': 6.5244, 'longitude': 3.3792},     # Nigeria (Esan region)
    'GWD': {'latitude': 13.4432, 'longitude': -16.7219},  # Gambia
    'JPT': {'latitude': 35.6762, 'longitude': 139.6503},  # Japan (Tokyo)
    'LWK': {'latitude': 0.0236, 'longitude': 34.7617},    # Kenya (Luhya region)
    'MSL': {'latitude': 8.4850, 'longitude': -13.2350},   # Sierra Leone
    'YRI': {'latitude': 7.3775, 'longitude': 3.9470}      # Nigeria (Yoruba region)
}

# Read your existing geo_data
geo_data = pd.read_csv('geographic_data.csv')

# Convert coordinates to numpy arrays right away
coordinates = np.array([[population_coords[pop]['latitude'], population_coords[pop]['longitude']] 
                       for pop in geo_data['population']], dtype=np.float64)

# Add as separate columns
geo_data['latitude'] = coordinates[:, 0]
geo_data['longitude'] = coordinates[:, 1]

# Verify the data types
print("Data types before saving:", geo_data.dtypes)

# Save with float specification
geo_data.to_csv('geo_data.csv', index=False, float_format='%.6f')

# Verify the file was saved correctly
test_read = pd.read_csv('geo_data.csv')
print("\nData types after reading saved file:", test_read.dtypes)
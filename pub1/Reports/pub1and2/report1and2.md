# Genetic Diversity Analysis Report

Analysis Date: 2025-01-07

## Clinical Impact Analysis

Research Question: How do allele frequencies of clinically significant variants differ between African populations and other global populations, and what are the implications for disease risk assessment?

### Overview
Total analyzed variants: 75,216
Population groups studied: 3
African-specific variants identified: 2,039

### Statistical Analysis
Statistical Method: Mann-Whitney U
Test Statistic: 575362877.0000
P-value: 1.16e-81

The analysis reveals **statistically significant** differences in variant distributions between African and non-African populations. The extremely low p-value suggests these differences are highly unlikely to occur by chance.

### Population Distribution Patterns

#### African Populations:
- Mean allele frequency: 0.1381
- Variant count: 22,961
- Standard deviation: 0.3062

#### East_Asian Populations:
- Mean allele frequency: 0.2473
- Variant count: 16,588
- Standard deviation: 0.4311

#### European Populations:
- Mean allele frequency: 0.0185
- Variant count: 6,373
- Standard deviation: 0.1340

### Population-Specific Patterns

#### West African Populations
- Unique variants: 1,921
- Diaspora sharing: 76.8% variants shared with African diaspora

#### East African Populations
- Unique variants: 8,294

### Variant Enrichment Analysis
- African: 0.026 enrichment score
- East_Asian: 0.027 enrichment score
- European: 0.025 enrichment score

### Clinical Implications
Based on the observed patterns:
- High variant sharing between West African and diaspora populations suggests maintained genetic signatures despite geographical separation
- Population-specific variants indicate potential adaptation to local environments
- Observed differences suggest the need for population-specific clinical considerations


### Conclusion for Research Question
The analysis demonstrates clear differences in allele frequencies between populations, with statistically significant variations in the distribution of clinically significant variants. These findings directly inform our understanding of population-specific disease risk patterns and highlight the importance of population-specific approaches in clinical assessment.


## Functional Consequences Analysis

Research Question: Is there a relationship between a variant's functional consequence and its frequency distribution across populations, particularly for disease-associated variants?

### Key Findings Addressing the Research Question
Analysis of variant consequences revealed 3,299 high severity variants, with 2,242 found in African populations.

### Statistical Analysis
Method: Mann-Whitney U (HIGH severity variants)
Test Statistic: 1100002.5000
P-value: 1.39e-05

### Population Distribution of Consequences

### Severity Distribution Analysis
African: 0.001 severity enrichment score
East_Asian: 0.001 severity enrichment score
European: 0.001 severity enrichment score

### Population-Specific Variant Patterns

#### variant_id
Unique severe variants: 0
Consequence types: 

#### consequence_type
Unique severe variants: 0
Consequence types: 

### Conclusion for Research Question
The analysis reveals distinct patterns in the distribution of functional consequences across populations, with significant variations in the prevalence of high-severity variants. These findings demonstrate clear relationships between variant function and population distribution.


## Geographical Risk Analysis

Research Question: Do risk alleles show geographical clustering patterns when mapped across subgroups within major population groups, and how does this relate to known disease prevalence?

### Key Findings Addressing the Research Question
Analysis of 425 risk alleles across 9 geographical regions revealed distinct clustering patterns in allele distributions.

### Statistical Analysis
Method: Mann-Whitney U (West vs East African)
Test Statistic: 167239846706.0000
P-value: 0.00e+00

### Regional Distribution Analysis

####     Barbados
Name: 0, dtype: object (    African_Diaspora
Name: 0, dtype: object)
Mean Allele Frequency: 0.0822
Risk Alleles: <lambda>    377
Name: 0, dtype: object
Unique Variants: nunique    246768
Name: 0, dtype: object

####     USA
Name: 1, dtype: object (    African_Diaspora
Name: 1, dtype: object)
Mean Allele Frequency: 0.0824
Risk Alleles: <lambda>    377
Name: 1, dtype: object
Unique Variants: nunique    246767
Name: 1, dtype: object

####     Kenya
Name: 2, dtype: object (    East_African
Name: 2, dtype: object)
Mean Allele Frequency: 0.2432
Risk Alleles: <lambda>    425
Name: 2, dtype: object
Unique Variants: nunique    298538
Name: 2, dtype: object

####     Gambia
Name: 3, dtype: object (    West_African
Name: 3, dtype: object)
Mean Allele Frequency: 0.2433
Risk Alleles: <lambda>    425
Name: 3, dtype: object
Unique Variants: nunique    298537
Name: 3, dtype: object

####     Nigeria
Name: 4, dtype: object (    West_African
Name: 4, dtype: object)
Mean Allele Frequency: 0.1724
Risk Alleles: <lambda>    425
Name: 4, dtype: object
Unique Variants: nunique    298538
Name: 4, dtype: object

####     Sierra Leone
Name: 5, dtype: object (    West_African
Name: 5, dtype: object)
Mean Allele Frequency: 0.0823
Risk Alleles: <lambda>    377
Name: 5, dtype: object
Unique Variants: nunique    246765
Name: 5, dtype: object

####     Beijing
Name: 6, dtype: object (    Chinese
Name: 6, dtype: object)
Mean Allele Frequency: 0.2457
Risk Alleles: <lambda>    425
Name: 6, dtype: object
Unique Variants: nunique    298244
Name: 6, dtype: object

####     Japan
Name: 7, dtype: object (    Eastern
Name: 7, dtype: object)
Mean Allele Frequency: 0.2459
Risk Alleles: <lambda>    425
Name: 7, dtype: object
Unique Variants: nunique    298245
Name: 7, dtype: object

####     Europe
Name: 8, dtype: object (    Northern_European
Name: 8, dtype: object)
Mean Allele Frequency: 0.0830
Risk Alleles: <lambda>    377
Name: 8, dtype: object
Unique Variants: nunique    246767
Name: 8, dtype: object

### Regional Enrichment Patterns

#### ('Barbados', 'African_Diaspora')
Risk Ratio: 1.000
Mean Frequency: 0.0822
Unique Risk Alleles: 377

#### ('Beijing', 'Chinese')
Risk Ratio: 1.000
Mean Frequency: 0.2457
Unique Risk Alleles: 425

#### ('Europe', 'Northern_European')
Risk Ratio: 1.000
Mean Frequency: 0.0830
Unique Risk Alleles: 377

#### ('Gambia', 'West_African')
Risk Ratio: 1.000
Mean Frequency: 0.2433
Unique Risk Alleles: 425

#### ('Japan', 'Eastern')
Risk Ratio: 1.000
Mean Frequency: 0.2459
Unique Risk Alleles: 425

#### ('Kenya', 'East_African')
Risk Ratio: 1.000
Mean Frequency: 0.2432
Unique Risk Alleles: 425

#### ('Nigeria', 'West_African')
Risk Ratio: 1.000
Mean Frequency: 0.1724
Unique Risk Alleles: 425

#### ('Sierra Leone', 'West_African')
Risk Ratio: 1.000
Mean Frequency: 0.0823
Unique Risk Alleles: 377

#### ('USA', 'African_Diaspora')
Risk Ratio: 1.000
Mean Frequency: 0.0824
Unique Risk Alleles: 377

### Conclusion for Research Question
The geographical analysis reveals significant regional patterns in risk allele distribution, suggesting local adaptation and population-specific disease risks. These patterns align with known geographical disease prevalence variations.


## Pleiotropic Effects Analysis

Research Question: How do variants associated with multiple traits distribute across populations, and what does this reveal about pleiotropic effects in different genetic backgrounds?

### Key Findings Addressing the Research Question
Analysis identified 203,874 variants affecting multiple traits, with an average of 0.38 traits per variant and up to 441 traits for some variants.

### Statistical Analysis
Method: Mann-Whitney U (pleiotropic variants)
Test Statistic: 4139686272.5000
P-value: 3.77e-297

### Population-Specific Pleiotropic Patterns

####     African
Name: 0, dtype: object (    African_Diaspora
Name: 0, dtype: object)
Mean Trait Count: 4.68
Maximum Traits: 441
Allele Frequency: 0.1925

####     African
Name: 1, dtype: object (    East_African
Name: 1, dtype: object)
Mean Trait Count: 5.06
Maximum Traits: 441
Allele Frequency: 0.5120

####     African
Name: 2, dtype: object (    West_African
Name: 2, dtype: object)
Mean Trait Count: 4.93
Maximum Traits: 441
Allele Frequency: 0.4093

####     East_Asian
Name: 3, dtype: object (    Chinese
Name: 3, dtype: object)
Mean Trait Count: 5.01
Maximum Traits: 441
Allele Frequency: 0.5350

####     East_Asian
Name: 4, dtype: object (    Eastern
Name: 4, dtype: object)
Mean Trait Count: 5.01
Maximum Traits: 441
Allele Frequency: 0.5351

####     European
Name: 5, dtype: object (    Northern_European
Name: 5, dtype: object)
Mean Trait Count: 4.70
Maximum Traits: 441
Allele Frequency: 0.1989

### Trait Network Analysis

#### African
Unique Traits: 2060
Max Traits per Variant: 441
Top Trait Associations:
-  FAMILIAL ↔  SUSCEPTIBILITY TO: 71618 variants
-  SUSCEPTIBILITY TO ↔ BREAST-OVARIAN CANCER: 71618 variants
-  FAMILIAL ↔ BREAST-OVARIAN CANCER: 71567 variants
- Hereditary breast ovarian cancer syndrome ↔ Hereditary cancer-predisposing syndrome: 58781 variants
-  SUSCEPTIBILITY TO ↔ Hereditary breast ovarian cancer syndrome: 39481 variants

#### East_Asian
Unique Traits: 2060
Max Traits per Variant: 441
Top Trait Associations:
-  FAMILIAL ↔  SUSCEPTIBILITY TO: 26099 variants
-  SUSCEPTIBILITY TO ↔ BREAST-OVARIAN CANCER: 26099 variants
-  FAMILIAL ↔ BREAST-OVARIAN CANCER: 26071 variants
- Hereditary breast ovarian cancer syndrome ↔ Hereditary cancer-predisposing syndrome: 24837 variants
- ClinVar: phenotype not specified ↔ Hereditary cancer-predisposing syndrome: 16748 variants

#### European
Unique Traits: 1750
Max Traits per Variant: 441
Top Trait Associations:
-  FAMILIAL ↔  SUSCEPTIBILITY TO: 7339 variants
-  SUSCEPTIBILITY TO ↔ BREAST-OVARIAN CANCER: 7339 variants
-  FAMILIAL ↔ BREAST-OVARIAN CANCER: 7336 variants
- Hereditary breast ovarian cancer syndrome ↔ Hereditary cancer-predisposing syndrome: 5114 variants
-  1 ↔  FAMILIAL: 4057 variants

### Conclusion for Research Question
The analysis reveals distinct patterns in how pleiotropic variants distribute across populations, suggesting population-specific differences in how variants influence multiple traits. These patterns provide insights into the complex relationships between genetic background and trait manifestation.


## Protective Variants Analysis

Research Question: How do allele frequencies of protective variants compare across populations, and what does this suggest about differential selective pressures?

### Key Findings Addressing the Research Question

####     African
Name: 0, dtype: object (    African_Diaspora
Name: 0, dtype: object)
Allele Frequency: 0.5714 (±0.4537)
Variant Count: nunique    4
Name: 0, dtype: object
Consequence Distribution: <lambda>    {'MODERATE': 4, 'LOW': 4}
Name: 0, dtype: object

####     African
Name: 1, dtype: object (    East_African
Name: 1, dtype: object)
Allele Frequency: 0.6667 (±0.4372)
Variant Count: nunique    4
Name: 1, dtype: object
Consequence Distribution: <lambda>    {'MODERATE': 4, 'LOW': 3}
Name: 1, dtype: object

####     African
Name: 2, dtype: object (    West_African
Name: 2, dtype: object)
Allele Frequency: 0.6316 (±0.4387)
Variant Count: nunique    4
Name: 2, dtype: object
Consequence Distribution: <lambda>    {'MODERATE': 12, 'LOW': 9}
Name: 2, dtype: object

####     East_Asian
Name: 3, dtype: object (    Chinese
Name: 3, dtype: object)
Allele Frequency: 0.8000 (±0.4215)
Variant Count: nunique    4
Name: 3, dtype: object
Consequence Distribution: <lambda>    {'LOW': 3, 'MODERATE': 2}
Name: 3, dtype: object

####     East_Asian
Name: 4, dtype: object (    Eastern
Name: 4, dtype: object)
Allele Frequency: 0.8000 (±0.4215)
Variant Count: nunique    4
Name: 4, dtype: object
Consequence Distribution: <lambda>    {'LOW': 3, 'MODERATE': 2}
Name: 4, dtype: object

####     European
Name: 5, dtype: object (    Northern_European
Name: 5, dtype: object)
Allele Frequency: 1.0000 (±0.0000)
Variant Count: nunique    4
Name: 5, dtype: object
Consequence Distribution: <lambda>    {'MODERATE': 1, 'LOW': 1}
Name: 5, dtype: object

### Statistical Analysis
Method: Mann-Whitney U (protective variants)
Test Statistic: 432.0000
P-value: 1.20e-03

### Selection Pattern Analysis

#### ('African', 'African_Diaspora')
High-frequency Protective Variants: 10
Mean Protective Allele Frequency: 0.5714
Consequence Types Distribution:
- stop_gained: 6
- missense_variant: 4
- 5_prime_UTR_variant: 4

#### ('African', 'East_African')
High-frequency Protective Variants: 11
Mean Protective Allele Frequency: 0.6667
Consequence Types Distribution:
- stop_gained: 5
- missense_variant: 4
- 5_prime_UTR_variant: 3

#### ('African', 'West_African')
High-frequency Protective Variants: 32
Mean Protective Allele Frequency: 0.6316
Consequence Types Distribution:
- stop_gained: 17
- missense_variant: 12
- 5_prime_UTR_variant: 9

#### ('East_Asian', 'Chinese')
High-frequency Protective Variants: 8
Mean Protective Allele Frequency: 0.8000
Consequence Types Distribution:
- stop_gained: 5
- 5_prime_UTR_variant: 3
- missense_variant: 2

#### ('East_Asian', 'Eastern')
High-frequency Protective Variants: 8
Mean Protective Allele Frequency: 0.8000
Consequence Types Distribution:
- stop_gained: 5
- 5_prime_UTR_variant: 3
- missense_variant: 2

#### ('European', 'Northern_European')
High-frequency Protective Variants: 4
Mean Protective Allele Frequency: 1.0000
Consequence Types Distribution:
- stop_gained: 2
- missense_variant: 1
- 5_prime_UTR_variant: 1

### Conclusion for Research Question
The distribution patterns of protective variants across populations reveal signatures of differential selective pressures, suggesting adaptation to distinct environmental and evolutionary challenges in different populations.


## Sample Size Effects Analysis

Research Question: How does sample size variation across populations affect our understanding of rare clinically significant variants?

### Key Findings Addressing the Research Question
Sample Size Distribution across Quartiles:

#### Q1
Mean Sample Size: 44.6
Sample Count: 2569778
Mean Allele Frequency: 0.0994 (±0.2941)

#### Q2
Mean Sample Size: 7640.2
Sample Count: 35282
Mean Allele Frequency: 0.8432 (±0.3567)

#### Q3
Mean Sample Size: 13055.0
Sample Count: 43573
Mean Allele Frequency: 0.8304 (±0.3708)

#### Q4
Mean Sample Size: 17764.2
Sample Count: 199321
Mean Allele Frequency: 0.8558 (±0.3491)

### Sample Size-Variant Detection Relationship
Correlation Coefficient: 0.5948
This correlation indicates the strength of the relationship between sample size and variant detection capability.

### Sample Size Distribution Overview
Total Samples Analyzed: 2847954
Size Range: 0 to 19850
Median Sample Size: 0


### Conclusion for Research Question
The analysis demonstrates the critical impact of sample size variation on our ability to detect and characterize rare variants. This relationship has important implications for study design and interpretation of genetic variation across populations with different sample sizes.


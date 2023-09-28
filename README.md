# Prioritization connectivity Colombia
Workflow to select complementary conservation areas based on different biodiversity attributes and connectivity metrics

## Data

* Ecoregions map
* Human footprint map
* RUNAP
* Any spatial layers used as features in the objective function. Examples: SDMs, Species Richness, Carbon index, water index
* Cost layer

## Codes

### Gap analysis

Code to calculate the percentage of the distribution area protected for each species and the total area that should be protected in the objective function

### Calculate Metrics Per PU

Code to generate planing units and calculate metrics per planning unit

### Ranking Planning Units Only

Code to rank planning units according to pareto function. The code performs the following actions:

1. Define objective function
2. Generate objective combinations
3. Apply objective function to all the scenario combinations
4. Ranks planning units 
5. Save planning unites

### Calcular conectividad

Calculate connectivity for different scenarios using Makurhini

## Authors and contributors

Andres Felipe Suarez Castro a.suarezcastro@griffith.edu.au

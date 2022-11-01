# Prioritization_connectivity_Colombia
Workflow to select complementary conservation areas based on different biodiversity attributes and connectivity metrics

## Data

Ecoregions map
Human footprint map
RUNAP
Any spatial layers used as features in the objective function. Examples: Species Richness, Carbon index, water index
Cost layer

## Codes

### calculate_metrics_per_pu

Code to generate planing units and calculate metrics per planning unit

### ranking_planning_units_only

Code to rank planning units according to pareto function. The code performs the following actions:

1. Define objective function
2. Generate objective combinations
3. Apply objective function to all the scenario combinations
4. Ranks planning units 
5. Save planning unites

### Calcular conectividad



# methylobacterium-arabidopsis-growth

This repository contains all data and scripts needed for reproducing all statistical analyses and figures of the article "Richness and composition of phyllosphere Methylobacterium communities cause variation in Arabidopsis thaliana growth".

## Notes on dataframes

### biom_data.csv

diversity : strains richness (number of strains) 

comm : community composition (each letter refer to a different strain; see below)

leaf_biom : leaf dry biomass (mg)

molds : presence of molds on plant and/or medium

shelf_x : position coordinates on shelves in growth chamber

shelf_y : position coordinates on shelves in growth chamber

### arabi_synthcomms.csv

A to L : each letter represents a Methylobacterium strain (see section below). In this dataframe, values are the proportion (relative abundance) of each corresponding strain in the community associated with a plant.

clades : phylogenetic clades of Methylobacterium strains (letters A, B, and D are unrelated to letters assigned to strains). These clades are based on this paper: 
Leducq, J. B., Sneddon, D., Santos, M., Condrain-Morel, D., Bourret, G., Martinez-Gomez, N. C., ... & Marx, C. J. (2022). Comprehensive phylogenomics of Methylobacterium reveals four evolutionary distinct groups and underappreciated phyllosphere diversity. Genome Biology and Evolution, 14(8), evac123.

### H2_all_models.csv

Parameters and statistical results for all models explored (n = 4096) in the context of our second hypothesis, based on the saturated model: dry leaf biomass ~ E-046 + J-078 + J-088 + J-059 + J-043 + J-067 + J-048 + J-076 + E-045 + J-092 + E-005 + J-068. Intercept indicates the mean biomass of samples when all strains included in the corresponding model are absent. ID: model number; +: parameters considered in the corresponding model; adj.R2: adjusted R2; df: degrees of freedom; LL: log-likelihood.

### H3_all_models.csv

Parameters and statistical results for all models explored (n = 167) in the context of our third hypothesis, based on the saturated model: dry leaf biomass ~ strain richness + J-067 + E-045 + J-092 + strain richness:J-067 + strain richness:E-045 + strain richness:J-092 + J-067:E-045 + J-067:J-092 + E-045:J-092 + strain richness:J-067:E-045 + strain richness:J-067:J-092 + strain richness:E-045:J-092 + J-067:E-045:J-092 + strain richness:J-067:E-045:J-092. Intercept indicates the mean biomass of one-strain samples that do not contain the strains included in the corresponding model. ID: model number; D2: polynomial strain richness; +: parameters considered in the corresponding model; adj.R2: adjusted R2; df: degrees of freedom; LL: log-likelihood.

## Methylobacterium strains

A: E-046

B: J-078

C: J-088

D: J-059

E: J-043

F: J-067

G: J-048

H: J-076

I: E-045

J: J-092

K: E-005

L: J-068

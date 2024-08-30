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

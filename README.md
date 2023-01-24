################### README
##### author : Gaëlle Munsch
##### date : 24/01/2023

These files are relative to the work presented in the preprint  Genome-wide association study of a semicontinuous trait: Illustration of the impact of the modeling strategy through the study of Neutrophil Extracellular Traps levels
Gaëlle Munsch, Carole Proust, Sylvie Labrouche-Colomer, Dylan Aïssi, Anne Boland, Pierre-Emmanuel Morange, Anne Roche, Luc de Chaisemartin, Annie Harroche, Robert Olaso, Jean-François Deleuze, Chloé James, Joseph Emmerich, David M Smadja, Hélène Jacqmin-Gadda, David-Alexandre Trégouët
medRxiv 2022.09.19.22279929; doi: https://doi.org/10.1101/2022.09.19.22279929

### Phenotypic_file.txt
This file contains the clinical variables used for the application. The variable "NET" has a semicontinuous distribution and represents the outcome variable. Other variables are covariates.

### RMSE.sh
This file contains the R code used for the three models of interest : Compound Poisson-Gamma, Negative Binomial and Tobit.
Bootstrap procedures were used to calculate the RMSE criteria with each model of interest.

### Simulation.sh
This file contains the R code used for simulations described in the paper. It includes the simulation of the genotypes, the model estimation (Compound Poisson-Gamma and Negative Binomial) and the number of significant tests.


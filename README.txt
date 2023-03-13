Description
==============

**Roilo, S., Paulus, A., Alarcón-Segura, V., Kock, L., Beckmann, M., Klein, N., Cord, A. A multidimensional approach is needed to better quantify land-use intensity in biodiversity models.
  Please acknowledge this paper when using the code or data outputs.


Roilo_2023_Land_use_intensity_metrics.R
--------------
This R script is used to run the analysis described in Roilo et al., to calculate and compare different land-use intensity (LUI) metrics covering the study region Mulde River Basin, Saxony, Germany.
The script is divided into the following sections:
 - grid preparation: describes the preparation of the square, hexagonal, and voronoi grids used in the study;
 - computation of land use variables: describes the calculation of each of the 10 LUI metrics;
 - mean LUI and SD across metrics: describes the calculation of mean LUI and standard deviation across the 10 metrics for each grid type;
 - make plots of the different grids & metrics: plots all grids and metrics.

Roilo_2023_LUI_virtual_species.R
--------------
This R script is used to create a virtual species with known land-use intensity (LUI)-species relationships to test how using different subsets of LUI metrics affects the outcomes of 
biodiversity models. For this last step, we used Generalised Additive Models to model the virtual species occurrence.
The script is divided into the following sections:
 - Create a virtual species: describes the generation of the virtual species based on selected LUI metrics;
 - Modelling virtual species occurrence: described the data preparation steps and the modelling framework using generalised additive models to explain the virtual species occurrence based on LUI metrics;
 - Ridgeplots of LUI values at the presence points of the virtual species: plots ridgeplots of the LUI values' distribution, as extracted at the presence points of the virtual species;
 - Distance measures across the density distributions: calculated the symmetrical Kullback-Leibler distance for each pair of metrics;
 - Correlation structures: plots correlation matrixes to explore the correlation structures of different LUI metrics.

Corresponding author: 
--------------
Stephanie Roilo **stephanie.roilo@tu-dresden.de**

Please get in touch if you have questions about the code or data.

For more information about the BESTMAP project, please see: www.bestmap.eu 

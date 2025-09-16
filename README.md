# Myelination-function-coupling
This repository provides core code and relevant toolboxes for data analysis in the article “Dual-axis myelination covariance drives the functional connectivity emergence during infancy”.

## This repository contains the following files:

### <0.Preprocessing>: 
The dHCP anatomical data surface-based preprocessing program 'dHCP_Term_anat.sh' includes registration, resampling and smoothing. dHCP fucntional data surface-based preprocessing program 'dHCP_Term_func.sh' includes volume-mapping, registration, resampling, smoothing and calculating functional connectome. After preprocessing, all data are ultimately alighed to the HCP-YA fs_LR space, and down-sampled to 5k_fs_LR mesh.

### <1.MFC_calculating>: 
Use 'MFC_calculating.m' to Calculate the vertex-level MFC/gMFC/sMFC for the demo dHCP subject (ID: sub-CC00056XX07 ses-10700), the demo data is provided in a subfolder 'data', and the results of demo data are stored in another subfolder 'result'.

### <2.Growth_effect_analysis>: 
Use 'Growth_effect_analysis.R' to investigate linear and nonlinear relationships between MFC and age, the required data is provided in a subfolder 'data', and the results are stored in another subfolder 'result'.

### <3.Distance_dependence_analysis>: 
Use 'Distance_dependence_analysis.m' to investigate the distance dependence of MFC.

### <4.Gene_analysis>: 
Use Brainspan Atlas data to investigate the gene association of MFC.

## The following analyses were carried out using open source packages:

### 1.Spin test: 
The spin test was conducted using an open Matlab code package (https://github.com/frantisekvasa/rotate_parcellation).

### 2.Generalized additive models (GAM): 
The generalized additive models were performed using the R package mgcv (https://cran.r-project.org/web/packages/mgcv/index.html).

### 3.Gradient calculation: 
The connectome gradient was computed using the Matlab BrainSpace package (https://www.mathworks.com/matlabcentral/fileexchange/86887-brainspace).

### 4.Gene Ontology Enrichment Analysis: 
Gene ontology enrichment analysis was conducted using Metascape (https://metascape.org/gp/index.html#/main/step1).

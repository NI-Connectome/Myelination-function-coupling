# Myelination-function-coupling
This repository provides core code and relevant toolboxes for data analysis in the article “Dual-axis myelination covariance drives the functional connectivity emergence during infancy”.

## This folder contains the following files:

### <0.Preprocessing>: 
The dHCP anatomical data surface-based preprocessing program 'dHCP_Term_anat.sh' includes registration, resampling and smoothing. dHCP fucntional data surface-based preprocessing program 'dHCP_Term_func.sh' includes volume-mapping, registration, resampling, smoothing and calculating functional connectome. After preprocessing, all data are ultimately alighed to the HCP-YA fs_LR space, and down-sampled to 5k_fs_LR mesh.

### <1.MFC_calculating>: 
Calculate the vertex-level MFC/gMFC/sMFC for 364 dHCP subjects ('MFC_calculating.m'), the information of 364 subjects is provided in the 'TermList_myelin.txt', and 'Label_7net_5k.mat' is Yeo's 7net label for surface-based analysis at 5k resolution.

### <2.Growth_effect_analysis>: 
Use GAM model to investigate linear and nonlinear relationships between MFC and age ('Growth_effect_analysis.R').

### <3.Distance_dependence_analysis>: 
Investigate the distance dependence of MFC ('Distance_dependence_analysis.m') by dividing vertices into multiple vertex-specific parts, then recalculating MFC for each part.

## The following analyses were carried out using open source packages:

### 1.Spin test: 
The spin test was conducted using an open Matlab code package (https://github.com/frantisekvasa/rotate_parcellation).

### 2.Generalized additive models (GAM): 
The generalized additive models were performed using the R package mgcv (https://cran.r-project.org/web/packages/mgcv/index.html).

### 3.Gradient calculation: 
The connectome gradient was computed using the Matlab BrainSpace package (https://www.mathworks.com/matlabcentral/fileexchange/86887-brainspace).

### 4.Gene Ontology Enrichment Analysis: 
Gene ontology enrichment analysis was conducted using Metascape (https://metascape.org/gp/index.html#/main/step1).

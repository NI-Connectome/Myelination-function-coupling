This folder contains the code required to calculate the growth effect of MFC, data to support the code, and the calculation results.

Specifically,

The program 'Growth_effect_analysis.R' uses the GAM model to capture the linear and nonlinear relationship between MFC and PMA and calculates the first derivative of the age smoothing term. The mean of the first derivative is used as the growth rate. 

The subfolder 'data' includes age (including GA, PMA and PNA of 364 subjects) and vertex-level MFC (the vertex-level MFC is too large to upload, please click the link below to download: https://pan.baidu.com/s/14Oi2W9z1itGbgNZqZexiBQ?pwd=jyzq).

The subfolder 'result' includes the growth rate measures as the mean of the first derivativecal (derivs_MFC.csv, derivs_gMFC.csv and derivs_sMFC.csv), and GAM model fitting effect measuread as p-value (Pvalue_MFC.csv, Pvalue_gMFC.csv and Pvalue_sMFC.cs).

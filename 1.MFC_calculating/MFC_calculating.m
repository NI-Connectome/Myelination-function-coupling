clc;
clear all;

% This code calculate the vertex-level MFC/gMFC/sMFC for the demo dHCP subject(ID:sub-CC00056XX07 ses-10700)
% The program includes calculating of FC, gMC, sMC, MFC, gMFC and sMFC, and visualization of brain map

your_path = './data'; % Please set your file storage path
load './data/Label_7net_5k.mat'  % Yeo's 7net label for surface-based analysis at 5k resolution
load './data/myelin_DemoSub.mat'  % myelinmap of demo_sub by downsampling to 5k resolution
load './data/myelin.mat'  % myelinmap of 364 subs by downsampling to 5k resolution
subj = 'CC00056XX07'; % the demo subject
sess = '10700'; % the demo session of the subject

Ind_notNuc = find(Label_7net_5k > 0);  % remove subcortical structures vertices
NumVertex = size(Ind_notNuc,1);     % number of cortical veretices: 8589
Ind_utri = find(triu(ones(NumVertex),1)); % matrix vectorization

%% load the calculated individual-specific FC for the demo subject
foldername = your_path;
filename = sprintf('sub-%s_ses-%s_hemi-LR_BOLD_correlation.hemi_5k.dconn.nii',subj,sess); % functional connectome
filepath = fullfile(foldername,filename);
sub_FC = ciftiopen(filepath);
sub_FC = sub_FC.cdata;

sub_FC = 0.5 * log((1+sub_FC) ./ (1-sub_FC));   % Fisher z-transform
sub_FC=sub_FC(Ind_notNuc,Ind_notNuc);   % 8589*8589 matrix

%% calculation of gMC (group-level myelination covariance)
myelin_z_Insub = zscore(myelin,0,1);  % z-score
gMC = corr(myelin_z_Insub');  % Pearson correlation inter-subject
gMC = 0.5 * log((1+gMC) ./ (1-gMC));  % Fisher z-transform

%% calculation of individual-specific sMC (subject-level myelination covariance) for the demo subject
sub_myelin = zscore(myelin_DemoSub);  % z-score
sub_MC = sub_myelin * sub_myelin';  % Pearson correlation intra-subject

%% calculation of individual-specific vertex-level gMFC for the demo subject
for m=1:length(Ind_notNuc)
    [b,~,r,~,stats]=regress(sub_FC(m,[1:m-1,m+1:end])',[ones(length(Ind_notNuc)-1,1),gMC(m,[1:m-1,m+1:end])']);  % FC~1+gMC
    MFC_gs_v_b(m,:) = b;
    sub_MFC_gs_v_re(m,:) = r;
    MFC_gs_v_stats(m,:) = stats;
    MFC_gs_v_R2(m,1) = stats(1);  
    MFC_gs_v_R2(m,1) = 1-(8588-1)/(8588-1-1)*(1-MFC_gs_v_R2(m,1)); % the adjusted coefficient of determination
end

%% calculation of individual-specific vertex-level sMFC for the demo subject
for m=1:length(Ind_notNuc)
    [b,~,r,~,stats]=regress(sub_FC(m,[1:m-1,m+1:end])',[ones(length(Ind_notNuc)-1,1),sub_MC(m,[1:m-1,m+1:end])']);
    MFC_ss_v_b(m,:) = b;
    sub_MFC_ss_v_re(m,:) = r;
    MFC_ss_v_stats(m,:) = stats;
    MFC_ss_v_R2(m,1) = stats(1);  
    MFC_ss_v_R2(m,1) = 1-(8588-1)/(8588-1-1)*(1-MFC_ss_v_R2(m,1)); % the adjusted coefficient of determination
end

%% calculation of individual-specific vertex-level MFC the demo subject
for m=1:length(Ind_notNuc)
    X_MCN = double([ones(length(Ind_notNuc)-1,1),zscore(gMC(m,[1:m-1,m+1:end])'),zscore(sub_MC(m,[1:m-1,m+1:end])')]);
    [b,~,r,~,stats]=regress(FC(m,[1:m-1,m+1:end])',X_MCN);
    MFC_gss_v_b(m,:) = b;
    sub_MFC_gss_v_re(m,:) = r;
    MFC_gss_v_stats(m,:) = stats;
    MFC_gss_v_R2(m,1) = stats(1); 
    MFC_gss_v_R2(m,1) = 1-(8588-1)/(8588-2-1)*(1-MFC_gss_v_R2(m,1)); % the adjusted coefficient of determination
end

%% draw brain map using workbench
mycifti = cifti_read('./data/sub-CC00056XX07_ses-10700_hemi-LR_sulc.hemi_5k.dtseries.nii');
cMFC_gs_v_R2 = zeros(9684,1);
cMFC_gs_v_R2(Ind_notNuc(1:8589,:),1) = MFC_gs_v_R2(1:8589,1);
mycifti.cdata = cMFC_gs_v_R2;  ciftisavereset(mycifti, './result/gMFC_DemoSub.dtseries.nii'); % gMFC

cMFC_ss_v_R2 = zeros(9684,1);
cMFC_ss_v_R2(Ind_notNuc(1:8589,:),1) = MFC_ss_v_R2(1:8589,1);
mycifti.cdata = cMFC_ss_v_R2;  ciftisavereset(mycifti, './result/sMFC_DemoSub.dtseries.nii'); % sMFC

cMFC_gss_v_R2 = zeros(9684,1);
cMFC_gss_v_R2(Ind_notNuc(1:8589,:),1) = MFC_gss_v_R2(1:8589,1);
mycifti.cdata = cMFC_gss_v_R2;  ciftisavereset(mycifti, './result/MFC_DemoSub.dtseries.nii'); % MFC

% Then, you can execute the 'wb_view' command in the system terminal to open the visualization interface, 
% import the surface file we provided (week-40_hemi-left_space-dhcpSym_dens-5k_midthickness.surf.gii and
% week-40_hemi-right_space-dhcpSym_dens-5k_midthickness.surf.gii) and the MFC data file obtained above 
% (MFC_DemoSub.dtseries.nii, gMFC_DemoSub.dtseries.nii, and sMFC_DemoSub.dtseries.nii), 
% which can display the MFC mind map at the vertex level. 
% The visualization effect is as shown in Figure 'DemoSub_brainmap_MFC.png'.
% You can adjust the colorbar according to your preferences.

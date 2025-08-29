clc;
clear all;

% this code calculate the vertex-level MFC/gMFC/sMFC for 364 dHCP subjects
% the information of 364 subjects is provided in the 'TermList_myelin.txt'

NumSub = 364; % number of subjects
load '/data/users/wliu/demo_dHCP_Analysis/Label_7net_5k.mat'  % Yeo's 7net label for surface-based analysis at 5k resolution
load '/data/users/wliu/demo_dHCP_Analysis/myelin.mat'  % myelinmap of 364 subs by downsampling to 5k resolution
Sublist = importdata('/data/users/wliu/demo_dHCP_Analysis/TermList_myelin.txt');  % 364 subs information
subj = string(Sublist.textdata);
sess = string(Sublist.data);

Ind_notNuc = find(Label_7net_5k > 0);  % remove subcortical structures vertices
NumVertex = size(Ind_notNuc,1);     % number of cortical veretices: 8589
Ind_utri = find(triu(ones(NumVertex),1)); % matrix vectorization

%% load rsfMRI data and calculate individual-specific FC
for i = 1:NumSub
    foldername = sprintf('/data/users/wliu/demo_dHCP_Analysis/sub-%s/ses-%s',subj(i),sess(i));
    filename = sprintf('sub-%s_ses-%s_hemi-LR_BOLD_correlation.hemi_5k.dconn.nii',subj(i),sess(i)); % functional connectome
    filepath = fullfile(foldername,filename);
    sub_FC = ciftiopen(filepath);
    sub_FC = sub_FC.cdata;

    sub_FC = 0.5 * log((1+sub_FC) ./ (1-sub_FC));   % Fisher z-transform
    sub_FC=sub_FC(Ind_notNuc,Ind_notNuc);   % 8589*8589 matrix
    subs_FC_Z_utri(:,i)=sub_FC_FisherZ(Ind_utri);   % matrix vectorization

    save(sprintf('/data/users/wliu/demo_dHCP_Analysis/sub-%s/ses-%s/sub_FC.mat',subj(i),sess(i)),'sub_FC');  % save the vectorized FC matrix
    fprintf('%dth subject is OK~\n',i);
end

%% calculation of gMC
myelin_z_Insub = zscore(myelin,0,1);  % z-score
gMC = corr(myelin_z_Insub');  % Pearson correlation inter-subject
gMC = 0.5 * log((1+gMC) ./ (1-gMC));  % Fisher z-transform

%% calculation of individual-specific sMC
for i=1:NumSub
    sub_myelin = zscore(myelin(:,i));  % z-score
    sub_MCN = sub_myelin * sub_myelin';  % Pearson correlation intra-subject
    sub_MCN = sub_MCN(Ind_utri);   % matrix vectorization
    
    save(sprintf('/data/users/wliu/demo_dHCP_Analysis/sub-%s/ses-%s/sub_MCN.mat',subj(i),sess(i)),'sub_MCN');  % save the vectorized sMC matrix
    fprintf('%dth subject is OK~\n',i);
end

%% calculation of individual-specific vertex-level gMFC
for k=1:NumSub
    FCname=sprintf('/data/users/wliu/demo_dHCP_Analysis/sub-%s/ses-%s/sub_FC.mat',subj(k),sess(k)); 
    load(FCname); % load individual-specific FC

    FC = zeros(length(Ind_notNuc),length(Ind_notNuc));
    FC(Ind_utri) = sub_FC;
    for i = 2:length(Ind_notNuc)
        for j = 1:i-1
            FC(i,j) = FC(j,i);   % 8589*8589 matrix
        end
    end

    for m=1:length(Ind_notNuc)
        [b,~,r,~,stats]=regress(FC(m,[1:m-1,m+1:end])',[ones(length(Ind_notNuc)-1,1),gMC(m,[1:m-1,m+1:end])']);  % FC~1+gMC
        MFC_gs_v_b(m,:,k) = b;
        sub_MFC_gs_v_re(m,:) = r;
        MFC_gs_v_stats(m,:,k) = stats;
        MFC_gs_v_R2(m,k) = stats(1);  
        MFC_ss_v_R2(m,k) = 1-(8588-1)/(8588-1-1)*(1-MFC_ss_v_R2(m,k))% the adjusted coefficient of determination
    end
    fprintf('%dth subject is OK~\n',k);
end

%% Vertex-level sMFC
for k=1:NumSub
    FCname=sprintf('/data/users/wliu/demo_dHCP_Analysis/sub-%s/ses-%s/sub_FC.mat',subj(k),sess(k));
    load(FCname);  % load individual-specific FC
    MCNname=sprintf('/data/users/wliu/demo_dHCP_Analysis/sub-%s/ses-%s/sub_MCN.mat',subj(k),sess(k));
    load(MCNname);  % load individual-specific sMC
 
    MCN = zeros(length(Ind_notNuc),length(Ind_notNuc));
    MCN(Ind_utri) = sub_MCN;
    for i = 2:length(Ind_notNuc)
        for j = 1:i-1
            MCN(i,j) = MCN(j,i);  % 8589*8589 sMC matrix
        end
    end

    FC = zeros(length(Ind_notNuc),length(Ind_notNuc));
    FC(Ind_utri) = sub_FC;
    for i = 2:length(Ind_notNuc)
        for j = 1:i-1
            FC(i,j) = FC(j,i);   % 8589*8589 FC matrix
        end 
    end

    for m=1:length(Ind_notNuc)
        [b,~,r,~,stats]=regress(FC(m,[1:m-1,m+1:end])',[ones(length(Ind_notNuc)-1,1),MCN(m,[1:m-1,m+1:end])']);
        MFC_ss_v_b(m,:,k) = b;
        sub_MFC_ss_v_re(m,:) = r;
        MFC_ss_v_stats(m,:,k) = stats;
        MFC_ss_v_R2(m,k) = stats(1);  
        MFC_ss_v_R2(m,k) = 1-(8588-1)/(8588-1-1)*(1-MFC_ss_v_R2(m,k))% the adjusted coefficient of determination
    end
    fprintf('%dth subject is OK~\n',k);
end

%% Vertex-level MFC
for k=1:NumSub
    FCname=sprintf('/data/users/wliu/demo_dHCP_Analysis/sub-%s/ses-%s/sub_FC.mat',subj(k),sess(k));
    load(FCname); % load individual-specific FC
    MCNname=sprintf('/data/users/wliu/demo_dHCP_Analysis/sub-%s/ses-%s/sub_MCN.mat',subj(k),sess(k));
    load(SCNname);  % load individual-specific sMC

    MCN = zeros(length(Ind_notNuc),length(Ind_notNuc));
    MCN(Ind_utri) = sub_MCN;
    for i = 2:length(Ind_notNuc)
        for j = 1:i-1
            MCN(i,j) = MCN(j,i);  % 8589*8589 sMC matrix
        end
    end

    FC = zeros(length(Ind_notNuc),length(Ind_notNuc));
    FC(Ind_utri) = sub_FC;
    for i = 2:length(Ind_notNuc)
        for j = 1:i-1
            FC(i,j) = FC(j,i);  % 8589*8589 FC matrix
        end
    end

    for m=1:length(Ind_notNuc)
        X_MCN = double([ones(length(Ind_notNuc)-1,1),zscore(gMC(m,[1:m-1,m+1:end])'),zscore(MCN(m,[1:m-1,m+1:end])')]);
        [b,~,r,~,stats]=regress(FC(m,[1:m-1,m+1:end])',X_MCN);
        MFC_gss_v_b(m,:,k) = b;
        sub_MFC_gss_v_re(m,:) = r;
        MFC_gss_v_stats(m,:,k) = stats;
        MFC_gss_v_R2(m,k) = stats(1); 
        MFC_gss_v_R2(m,k) = 1-(8588-1)/(8588-2-1)*(1-MFC_gss_v_R2(m,k))% the adjusted coefficient of determination
    end
    fprintf('%dth subject is OK~\n',k);
end


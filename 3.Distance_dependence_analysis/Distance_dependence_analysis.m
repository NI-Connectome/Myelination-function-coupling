clc;
clear all;

% this code calculate the distance dependence of MFC
% the geodesic distance betweeen vertices is measured on the mnid-thickness surface

subj = 'CC00056XX07'; % the demo subject
sess = '10700'; % the demo session of the subject

%%%% load geodesic distance matrix on the mid-thickness surafce 
load('./data/distance.mat','distance_left')  % geodesic distance matrix for left sphere
load('./data/distance.mat','distance_right')  % geodesic distance matrix for right sphere

%% divide vertices into 50 vertex-specific parts of equal size based on vertex-to-vertex distance
%%%%% vertex-specific 50 parts for left sphere
for i=1:50 
    Ind_sort(i,1)=floor(2+85.8*(i-1));
    Ind_sort(i,2)=floor(2+85.8*i)-1;  % 42.9
end
for i=1:size(distance_left,1)
    distance_node = distance_left(i,:);
    [~,sortedInd] = sort(distance_node,'descend');
    for j=1:50
        subnode_left{i,j} = sortedInd(Ind_sort(j,1):Ind_sort(j,2));  % vertex-specific 50 parts
    end
end
%%%%% vertex-specific 50 parts for right sphere
for i=1:50
    Ind_sort(i,1)=floor(2+85.94*(i-1));
    Ind_sort(i,2)=floor(2+85.94*i)-1;    % 42.97
end
for i=1:size(distance_right,1)
    distance_node = distance_right(i,:);
    [~,sortedInd] = sort(distance_node,'descend');
    for j=1:50
        subnode_right{i,j} = sortedInd(Ind_sort(j,1):Ind_sort(j,2));  % vertex-specific 50 parts
    end
end

%% calculate MFC of vertex-specific 50 parts
for i=1:size(distance_left,1)
    for j=1:50
    subnode_left_now = subnode_left{i,j};  % vertex-specific part
    MCNg_left = gMC(i,subnode_left_now)';
    MCNs_left = sub_MC(i,subnode_left_now)';
    FCs_left = sub_FC(i,subnode_left_now)';
    X_MCN = double([ones(size(MCNg_left,1),1),zscore(MCNg_left),zscore(MCNs_left)]);
    [~,~,~,~,stats] = regress(FCs_left,X_MCN);   % MFC of each part
    MFC_gss_left(i,j) = stats(1);
    MFC_gss_left(i,j) = 1-(8588-1)/(8588-2-1)*(1-MFC_gss_left(i,j)); % the adjusted coefficient of determination
    end
end

for i=1:size(distance_right,1)
    for j=1:50
    subnode_right_now = subnode_right{i,j};   % vertex-specific part
    MCNg_right = gMC(i+4291,subnode_right_now+4291)';
    MCNs_right = sub_MC(i+4291,subnode_right_now+4291)';
    FCs_right = sub_FC(i+4291,subnode_right_now+4291)';
    X_MCN = double([ones(size(MCNg_right,1),1),zscore(MCNg_right),zscore(MCNs_right)]);
    [~,~,~,~,stats] = regress(FCs_right,X_MCN);    % MFC of each part
    MFC_gss_right(i,j) = stats(1);
    MFC_gss_right(i,j) = 1-(8588-1)/(8588-2-1)*(1-MFC_gss_right(i,j)); % the adjusted coefficient of determination
    end
end


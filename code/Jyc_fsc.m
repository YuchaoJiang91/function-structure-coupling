%% compute nodal function-structure coupling (fsc) using Spearman rank correlation between FC and SC of each node with all other nodes
clear;clc

%% data for test 
%cd('F:\cheng_fsc\code')
%load('test_data.mat');
%sc = test_sc_count; fc = test_fc;
%clearvars -except sc fc;

%% data for ukb
cd('F:\cheng_fsc\Paired')
load('aparc_S1_FC_pearson_z_GSR.mat');
fc = FC; 
load('aparc_S1_SC_count.mat');
sc = SC_count;
clearvars -except sc fc;

%% FSC calculation
% set weak structural connection to zero
sc_threshold = 3;   
sc(sc<sc_threshold) = 0;

% calculate nodal fsc for each subject
[n_roi,~,n_subj] = size(fc);
fsc = zeros(n_subj,n_roi);
sc_edge_number = zeros(n_subj,n_roi);
for i = 1:n_subj
    ii = i./1000;
    if ii == fix(ii)
        disp(['n_subj = ', num2str(i)])
    end
    subj_all_sc = sc(:,:,i);
    subj_all_fc = fc(:,:,i); 
    subj_all_sc = triu(subj_all_sc,1) + tril(subj_all_sc,-1);  % set diagonal to 0
    subj_all_fc = triu(subj_all_fc,1) + tril(subj_all_fc,-1); 
    for k = 1:n_roi
        subj_node_sc = subj_all_sc(:,k);
        subj_node_fc = subj_all_fc(:,k);
        subj_node_sc_nonzero = subj_node_sc(subj_node_sc>0);
        subj_node_fc_nonzero = subj_node_fc(subj_node_sc>0);        
        sc_edge_number(i,k) = length(subj_node_sc_nonzero);        
        if sc_edge_number(i,k) ~= 0
            [r,p] = corr(subj_node_sc_nonzero,subj_node_fc_nonzero,'type','spearman');
            fsc(i,k) = r;
        end
    end
end
fsc_global = mean(fsc,2);
clearvars -except  sc fc fsc fsc_global sc_edge_number 
save('ukb_fsc.mat','fsc','fsc_global','sc_edge_number');
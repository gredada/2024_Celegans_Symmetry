clear;clc;

%% Mirror Symmetry
clear;clc;
cwd = pwd;
load('data/celegans_connectome.mat');

[Herm_Mirror_Symmetry,Herm_Mirror_Symmetry_null]=asymtool_mirror_symmetry(Herm_Adj,Herm_cont_ind);
[Male_Mirror_Symmetry,Male_Mirror_Symmetry_null]=asymtool_mirror_symmetry(Male_Adj,Male_cont_ind);

mkdir('results');
cd results
save mirror_symmetry_index *Mirror*
cd(cwd);


%% Mirror Symmetry display
iter_num=length(Herm_Mirror_Symmetry_null);
data_plot = [Herm_Mirror_Symmetry,mean(Herm_Mirror_Symmetry_null);...
             Male_Mirror_Symmetry,mean(Male_Mirror_Symmetry_null)];
bar(data_plot);
xticks(1:2);
xticklabels({'Hermaphrodite','Male'});
title('Mirror-symmetry index');  
hold on;
ngroups = size(data_plot, 1); nbars = size(data_plot, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
if ~floor(iter_num*0.025)
    iter_num=1/0.025;
end

% % % get 95% CI % % %
temp_max_herm_null = maxk((Herm_Mirror_Symmetry_null),floor(iter_num*0.025));
temp_min_herm_null = mink((Herm_Mirror_Symmetry_null),floor(iter_num*0.025));
temp_max_male_null = maxk((Male_Mirror_Symmetry_null),floor(iter_num*0.025));
temp_min_male_null = mink((Male_Mirror_Symmetry_null),floor(iter_num*0.025));

err_min = -[temp_min_herm_null(end)-mean(Herm_Mirror_Symmetry_null);...
           temp_min_male_null(end)-mean(Male_Mirror_Symmetry_null)];
err_max = [temp_max_herm_null(end)-mean(Herm_Mirror_Symmetry_null);...
           temp_max_male_null(end)-mean(Male_Mirror_Symmetry_null)];
x = (1:ngroups) - groupwidth/2 + (2*2-1) * groupwidth / (2*nbars);
errorbar(x, data_plot(:,2), err_min(:,1),err_max(:,1), '.');
legend({'C. elegans data','Null model','error bar'},'Location','northeast');
ylim([0,1])
hold off

%% pair-robustness
clear;clc;
cwd = pwd;
load('data/celegans_connectome.mat');

[Herm_Pair_Redundancy_bilateral,Herm_Pair_Redundancy_all] = asymtool_pair_redundancy_index(Herm_Adj,Herm_LRU);
[Male_Pair_Redundancy_bilateral,Male_Pair_Redundancy_all] = asymtool_pair_redundancy_index(Male_Adj,Male_LRU);

cd results
save pair_redundancy_index *Pair_Redundancy*
cd(cwd);

% p_value
p_value_herm=custom_shuffle_test(Herm_Pair_Redundancy_all,Herm_Pair_Redundancy_bilateral,100000,1);
p_value_male=custom_shuffle_test(Male_Pair_Redundancy_all,Male_Pair_Redundancy_bilateral,100000,1);
disp(p_value_herm)
disp(p_value_male)
%% null model for pair-robustness
clear;clc;
cwd = pwd;
load('data/celegans_connectome.mat');

null_n=1000;
Herm_pair_n = sum(Herm_LRU(:,1));
Male_pair_n = sum(Male_LRU(:,1));

Herm_Pair_Redundancy_bilateral_null = zeros(Herm_pair_n,null_n);
Male_Pair_Redundancy_bilateral_null = zeros(Male_pair_n,null_n);

for i=1:null_n
    disp(['generating random network ',num2str(i),])
    
    null_Herm_Adj=randmio_dir(Herm_Adj,1000);
    null_Male_Adj=randmio_dir(Male_Adj,1000);
    
    [Herm_Pair_Redundancy_bilateral_null(:,i),~] ...
        = asymtool_pair_redundancy_index(null_Herm_Adj,Herm_LRU,true,false);
    [Male_Pair_Redundancy_bilateral_null(:,i),~] ...
        = asymtool_pair_redundancy_index(null_Male_Adj,Male_LRU,true,false);
end

cd results
save pair_redundancy_index_null *Pair_Redundancy*
cd(cwd);

%% display for pair-robustness
% figure; plotSpread({Herm_Pair_Redundancy_bilateral,Herm_Pair_Redundancy_all,Male_Pair_Redundancy_bilateral,Male_Pair_Redundancy_all},'showMM',5,'distributionColor','k')

%% connectivity similarity
clear; clc;
cwd = pwd;
load('data/celegans_connectome.mat');

% jaccard index (1-step)
[Herm_Jaccard_Index_1step_in_out_bilateral, Herm_Jaccard_Index_1step_in_out_all, ...
    Herm_Jaccard_Index_1step_in_bilateral, Herm_Jaccard_Index_1step_in_all] ...
    = asymtool_jaccard_1step(Herm_Adj,Herm_LRU);
[Male_Jaccard_Index_1step_in_out_bilateral, Male_Jaccard_Index_1step_in_out_all, ...
    Male_Jaccard_Index_1step_in_bilateral, Male_Jaccard_Index_1step_in_all] ...
    = asymtool_jaccard_1step(Male_Adj,Male_LRU);

% jaccard index (2-step)
[Herm_Jaccard_Index_2step_in_out_bilateral, Herm_Jaccard_Index_2step_in_out_all, ...
    Herm_Jaccard_Index_2step_in_bilateral, Herm_Jaccard_Index_2step_in_all] ...
    = asymtool_jaccard_2step(Herm_Adj,Herm_LRU);
[Male_Jaccard_Index_2step_in_out_bilateral, Male_Jaccard_Index_2step_in_out_all, ...
    Male_Jaccard_Index_2step_in_bilateral, Male_Jaccard_Index_2step_in_all] ...
    = asymtool_jaccard_2step(Male_Adj,Male_LRU);

cd results
save jaccard_index *Jaccard_Index*
cd(cwd);

% p_value
p_value_herm=custom_shuffle_test(Herm_Jaccard_Index_1step_in_out_all,Herm_Jaccard_Index_1step_in_out_bilateral,100000,1);
p_value_male=custom_shuffle_test(Male_Jaccard_Index_2step_in_out_all,Male_Jaccard_Index_2step_in_out_bilateral,100000,1);
disp(p_value_herm)
disp(p_value_male)
%% display for connectivity similarity

% figure; plotSpread({Herm_Jaccard_Index_1step_in_out_bilateral,Herm_Jaccard_Index_1step_in_out_all},'showMM',5,'distributionColor','k')
% figure; plotSpread({Herm_Jaccard_Index_1step_in_bilateral,Herm_Jaccard_Index_1step_in_all},'showMM',5,'distributionColor','k')
% figure; plotSpread({Male_Jaccard_Index_1step_in_out_bilateral,Male_Jaccard_Index_1step_in_out_all},'showMM',5,'distributionColor','k')
% figure; plotSpread({Male_Jaccard_Index_1step_in_bilateral,Male_Jaccard_Index_1step_in_all},'showMM',5,'distributionColor','k')
% figure; plotSpread({Herm_Jaccard_Index_2step_in_out_bilateral,Herm_Jaccard_Index_2step_in_out_all},'showMM',5,'distributionColor','k')
% figure; plotSpread({Herm_Jaccard_Index_2step_in_bilateral,Herm_Jaccard_Index_2step_in_all},'showMM',5,'distributionColor','k')
% figure; plotSpread({Male_Jaccard_Index_2step_in_out_bilateral,Male_Jaccard_Index_2step_in_out_all},'showMM',5,'distributionColor','k')
% figure; plotSpread({Male_Jaccard_Index_2step_in_bilateral,Male_Jaccard_Index_2step_in_all},'showMM',5,'distributionColor','k')
%% null model for connectivity similarity
clear;clc;
cwd = pwd;
load('data/celegans_connectome.mat');

null_n=1000;
Herm_pair_n = sum(Herm_LRU(:,1));
Male_pair_n = sum(Male_LRU(:,1));

Herm_Jaccard_Index_1step_in_out_bilateral_null  = zeros(Herm_pair_n,null_n);
Herm_Jaccard_Index_1step_in_bilateral_null      = zeros(Herm_pair_n,null_n);
Male_Jaccard_Index_1step_in_out_bilateral_null  = zeros(Male_pair_n,null_n);
Male_Jaccard_Index_1step_in_bilateral_null      = zeros(Male_pair_n,null_n);

Herm_Jaccard_Index_2step_in_out_bilateral_null  = zeros(Herm_pair_n,null_n);
Herm_Jaccard_Index_2step_in_bilateral_null      = zeros(Herm_pair_n,null_n);
Male_Jaccard_Index_2step_in_out_bilateral_null  = zeros(Male_pair_n,null_n);
Male_Jaccard_Index_2step_in_bilateral_null      = zeros(Male_pair_n,null_n);

for i=1:null_n
    disp(['generating random network ',num2str(i),])
    
    null_Herm_Adj=randmio_dir(Herm_Adj,1000);
    null_Male_Adj=randmio_dir(Male_Adj,1000);
    
    % jaccard index (1-step)
    [Herm_Jaccard_Index_1step_in_out_bilateral_null(:,i), ~, ...
        Herm_Jaccard_Index_1step_in_bilateral_null(:,i), ~] ...
        = asymtool_jaccard_1step(null_Herm_Adj,Herm_LRU);
    [Male_Jaccard_Index_1step_in_out_bilateral_null(:,i), ~, ...
        Male_Jaccard_Index_1step_in_bilateral_null(:,i), ~] ...
        = asymtool_jaccard_1step(null_Male_Adj,Male_LRU);
    
    % jaccard index (2-step)
    [Herm_Jaccard_Index_2step_in_out_bilateral_null(:,i), ~, ...
        Herm_Jaccard_Index_2step_in_bilateral_null(:,i), ~] ...
        = asymtool_jaccard_2step(null_Herm_Adj,Herm_LRU);
    [Male_Jaccard_Index_2step_in_out_bilateral_null(:,i), ~, ...
        Male_Jaccard_Index_2step_in_bilateral_null(:,i), ~] ...
        = asymtool_jaccard_2step(null_Male_Adj,Male_LRU);
end

cd results
save jaccard_index_null *Jaccard_Index*
cd(cwd);


%% Symmetry/Asymmetry breaking for connectivity similarity
clear;clc;
cwd = pwd;
load('data/celegans_connectome.mat');

null_n=1000;
Herm_pair_n = sum(Herm_LRU(:,1));
Male_pair_n = sum(Male_LRU(:,1));

Herm_Jaccard_Index_1step_in_out_bilateral_asymlist  = zeros(Herm_pair_n,null_n);
Herm_Jaccard_Index_1step_in_bilateral_asymlist      = zeros(Herm_pair_n,null_n);
Male_Jaccard_Index_1step_in_out_bilateral_asymlist  = zeros(Male_pair_n,null_n);
Male_Jaccard_Index_1step_in_bilateral_asymlist      = zeros(Male_pair_n,null_n);
Herm_Jaccard_Index_2step_in_out_bilateral_asymlist  = zeros(Herm_pair_n,null_n);
Herm_Jaccard_Index_2step_in_bilateral_asymlist      = zeros(Herm_pair_n,null_n);
Male_Jaccard_Index_2step_in_out_bilateral_asymlist  = zeros(Male_pair_n,null_n);
Male_Jaccard_Index_2step_in_bilateral_asymlist      = zeros(Male_pair_n,null_n);

Herm_Jaccard_Index_1step_in_out_bilateral_symlist  = zeros(Herm_pair_n,null_n);
Herm_Jaccard_Index_1step_in_bilateral_symlist      = zeros(Herm_pair_n,null_n);
Male_Jaccard_Index_1step_in_out_bilateral_symlist  = zeros(Male_pair_n,null_n);
Male_Jaccard_Index_1step_in_bilateral_symlist      = zeros(Male_pair_n,null_n);
Herm_Jaccard_Index_2step_in_out_bilateral_symlist  = zeros(Herm_pair_n,null_n);
Herm_Jaccard_Index_2step_in_bilateral_symlist      = zeros(Herm_pair_n,null_n);
Male_Jaccard_Index_2step_in_out_bilateral_symlist  = zeros(Male_pair_n,null_n);
Male_Jaccard_Index_2step_in_bilateral_symlist      = zeros(Male_pair_n,null_n);

for p_asym=0:0.1:1
    for i=1:1000
        temp_A=randmio_dir_ratio(Herm_Adj, Herm_cont_ind, p_asym, 0, 1000);
        temp_B=randmio_dir_ratio(Male_Adj, Male_cont_ind, p_asym, 0, 1000);
        % jaccard index (1-step)
        [Herm_Jaccard_Index_1step_in_out_bilateral_asymlist(:,i), ~, ...
            Herm_Jaccard_Index_1step_in_bilateral_asymlist(:,i), ~] ...
            = asymtool_jaccard_1step(temp_A,Herm_LRU);
        [Male_Jaccard_Index_1step_in_out_bilateral_asymlist(:,i), ~, ...
            Male_Jaccard_Index_1step_in_bilateral_asymlist(:,i), ~] ...
            = asymtool_jaccard_1step(temp_B,Male_LRU);
        
        % jaccard index (2-step)
        [Herm_Jaccard_Index_2step_in_out_bilateral_asymlist(:,i), ~, ...
            Herm_Jaccard_Index_2step_in_bilateral_asymlist(:,i), ~] ...
            = asymtool_jaccard_2step(temp_A,Herm_LRU);
        [Male_Jaccard_Index_2step_in_out_bilateral_asymlist(:,i), ~, ...
            Male_Jaccard_Index_2step_in_bilateral_asymlist(:,i), ~] ...
            = asymtool_jaccard_2step(temp_B,Male_LRU);
    end
end
for p_sym=0:0.1:1
    for i=1:1000
        temp_A=randmio_dir_ratio(Herm_Adj, Herm_cont_ind, 0, p_sym, 1000);
        temp_B=randmio_dir_ratio(Male_Adj, Male_cont_ind, 0, p_sym, 1000);
        % jaccard index (1-step)
        [Herm_Jaccard_Index_1step_in_out_bilateral_symlist(:,i), ~, ...
            Herm_Jaccard_Index_1step_in_bilateral_symlist(:,i), ~] ...
            = asymtool_jaccard_1step(temp_A,Herm_LRU);
        [Male_Jaccard_Index_1step_in_out_bilateral_symlist(:,i), ~, ...
            Male_Jaccard_Index_1step_in_bilateral_symlist(:,i), ~] ...
            = asymtool_jaccard_1step(temp_B,Male_LRU);
        
        % jaccard index (2-step)
        [Herm_Jaccard_Index_2step_in_out_bilateral_symlist(:,i), ~, ...
            Herm_Jaccard_Index_2step_in_bilateral_symlist(:,i), ~] ...
            = asymtool_jaccard_2step(temp_A,Herm_LRU);
        [Male_Jaccard_Index_2step_in_out_bilateral_symlist(:,i), ~, ...
            Male_Jaccard_Index_2step_in_bilateral_symlist(:,i), ~] ...
            = asymtool_jaccard_2step(temp_B,Male_LRU);
    end
end

cd results
save jaccard_index_break *Jaccard_Index*
cd(cwd);

%     data = jaccard_index_1step_in_Herm_sym_list';
%     
%     temp=mink(data,2);
%     errmin=temp(end,:)-mean(data);
%     temp=maxk(data,2);
%     errmax=temp(end,:)-mean(data);
%     errorbar(symm_breaking_ratio_list, mean(data), errmin,errmax, '-k.');
%% motif fingerprint difference
clear; clc;
cwd = pwd;
load('data/celegans_connectome.mat');

% Motif-fingerprint difference (3 nodes)
[Herm_MFdiff_3node_struc_bilateral,Herm_MFdiff_3node_struc_all, ...
    Herm_MFdiff_3node_func_bilateral,Herm_MFdiff_3node_func_all] ...
    = asymtool_motif_fingerprint_difference(Herm_Adj,Herm_LRU,3);
[Male_MFdiff_3node_struc_bilateral,Male_MFdiff_3node_struc_all, ...
    Male_MFdiff_3node_func_bilateral,Male_MFdiff_3node_func_all] ...
    = asymtool_motif_fingerprint_difference(Male_Adj,Male_LRU,3);

% Motif-fingerprint difference (4 nodes)
[Herm_MFdiff_4node_struc_bilateral,Herm_MFdiff_4node_struc_all, ...
    Herm_MFdiff_4node_func_bilateral,Herm_MFdiff_4node_func_all] ...
    = asymtool_motif_fingerprint_difference(Herm_Adj,Herm_LRU,4);
[Male_MFdiff_4node_struc_bilateral,Male_MFdiff_4node_struc_all, ...
    Male_MFdiff_4node_func_bilateral,Male_MFdiff_4node_func_all] ...
    = asymtool_motif_fingerprint_difference(Male_Adj,Male_LRU,4);

cd results
save motif_fingerprint_difference *MFdiff*
cd(cwd);

% p_value
p_value_herm=custom_shuffle_test(Herm_MFdiff_3node_struc_all,Herm_MFdiff_3node_struc_bilateral,100000,0);
p_value_male=custom_shuffle_test(Male_MFdiff_4node_func_all,Male_MFdiff_4node_func_bilateral,100000,0);
disp(p_value_herm)
disp(p_value_male)
%% null model for motif fingerprint difference
clear;clc;
cwd = pwd;
load('data/celegans_connectome.mat');

null_n=1000;
Herm_pair_n = sum(Herm_LRU(:,1));
Male_pair_n = sum(Male_LRU(:,1));

Herm_MFdiff_3node_struc_bilateral_null = zeros(Herm_pair_n,null_n);
Herm_MFdiff_3node_func_bilateral_null = zeros(Herm_pair_n,null_n);
Male_MFdiff_3node_struc_bilateral_null = zeros(Male_pair_n,null_n);
Male_MFdiff_3node_func_bilateral_null = zeros(Male_pair_n,null_n);

Herm_MFdiff_4node_struc_bilateral_null = zeros(Herm_pair_n,null_n);
Herm_MFdiff_4node_func_bilateral_null = zeros(Herm_pair_n,null_n);
Male_MFdiff_4node_struc_bilateral_null = zeros(Male_pair_n,null_n);
Male_MFdiff_4node_func_bilateral_null = zeros(Male_pair_n,null_n);

for i=1:null_n
    disp(['generating random network ',num2str(i),])
    
    null_Herm_Adj=randmio_dir(Herm_Adj,1000);
    null_Male_Adj=randmio_dir(Male_Adj,1000);
    
    [Herm_MFdiff_3node_struc_bilateral_null(:,i), ~, Herm_MFdiff_3node_func_bilateral_null(:,i), ~] ...
        = asymtool_motif_fingerprint_difference(null_Herm_Adj,Herm_LRU,3,true,false,true,true);
    [Herm_MFdiff_4node_struc_bilateral_null(:,i), ~, Herm_MFdiff_4node_func_bilateral_null(:,i), ~] ...
        = asymtool_motif_fingerprint_difference(null_Herm_Adj,Herm_LRU,4,true,false,true,true);
    
    [Male_MFdiff_3node_struc_bilateral_null(:,i), ~, Male_MFdiff_3node_func_bilateral_null(:,i), ~] ...
        = asymtool_motif_fingerprint_difference(null_Male_Adj,Male_LRU,3,true,false,true,true);
    [Male_MFdiff_4node_struc_bilateral_null(:,i), ~, Male_MFdiff_4node_func_bilateral_null(:,i), ~] ...
        = asymtool_motif_fingerprint_difference(null_Male_Adj,Male_LRU,4,true,false,true,true);    

end

cd results
save motif_fingerprint_difference_null *MFdiff*
cd(cwd);
%% display for motif fingerprint difference
clear; clc; close all;
cd results
load motif_fingerprint_difference_null
load motif_fingerprint_difference
cd ..
for i_txt=[{'3node_func'},{'3node_struc'},{'4node_func'},{'4node_struc'}]

eval([ 'temp_herm_all               = Herm_MFdiff_',i_txt{:},'_all;'            ])
eval([ 'temp_herm_bilateral         = Herm_MFdiff_',i_txt{:},'_bilateral;'      ])
eval([ 'temp_herm_bilateral_null    = Herm_MFdiff_',i_txt{:},'_bilateral_null;' ])
eval([ 'temp_male_all               = Male_MFdiff_',i_txt{:},'_all;'            ])
eval([ 'temp_male_bilateral         = Male_MFdiff_',i_txt{:},'_bilateral;'      ])
eval([ 'temp_male_bilateral_null    = Male_MFdiff_',i_txt{:},'_bilateral_null;' ])

figure;
data_plot = [mean(temp_herm_all),mean(temp_herm_bilateral),mean(mean(temp_herm_bilateral_null));...
    mean(temp_male_all),mean(temp_male_bilateral),mean(mean(temp_male_bilateral_null))];
b=bar(data_plot);
xticks([1,2]);
xticklabels({'hermaphrodite','male'});
title(['MFdiff for ',i_txt{:}(1:5),' ',i_txt{:}(7:end)]);   %
% error bars for null
hold on;
null_n=size(temp_herm_bilateral_null,2);
ngroups = size(data_plot, 1); nbars = size(data_plot, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
if ~floor(null_n*0.025)
    null_n=1/0.025;
end

% % % get 95% CI % % %
temp_max_struc_null = maxk(mean(temp_herm_bilateral_null),floor(null_n*0.025));
temp_min_struc_null = mink(mean(temp_herm_bilateral_null),floor(null_n*0.025));
temp_max_func_null = maxk(mean(temp_male_bilateral_null),floor(null_n*0.025));
temp_min_func_null = mink(mean(temp_male_bilateral_null),floor(null_n*0.025));

err_min = -[temp_min_struc_null(end)-mean(mean(temp_herm_bilateral_null));temp_min_func_null(end)-mean(mean(temp_male_bilateral_null))];
err_max = [temp_max_struc_null(end)-mean(mean(temp_herm_bilateral_null));temp_max_func_null(end)-mean(mean(temp_male_bilateral_null))];
x = (1:ngroups) - groupwidth/2 + (2*3-1) * groupwidth / (2*nbars);
e=errorbar(x, data_plot(:,3), err_min(:,1),err_max(:,1), 'r.');
b(3).FaceColor=[1 1 1];
b(2).FaceColor=[.5 .5 .5];
b(1).FaceColor=[0 0 0];
legend({'All','Bilateral','null'},'Location','northeast');
% ylim([0,1])
hold off
end
%% PCA group analysis
clear; clc;
close all

% data
load('data/celegans_connectome.mat');

load('results/pair_redundancy_index.mat');
load('results/jaccard_index.mat');
load('results/motif_fingerprint_difference.mat');

herm_pr_data = Herm_Pair_Redundancy_bilateral;
male_pr_data = Male_Pair_Redundancy_bilateral;

herm_cs_data = [Herm_Jaccard_Index_1step_in_out_bilateral, Herm_Jaccard_Index_1step_in_bilateral, ...
    Herm_Jaccard_Index_2step_in_out_bilateral, Herm_Jaccard_Index_2step_in_bilateral];
male_cs_data = [Male_Jaccard_Index_1step_in_out_bilateral, Male_Jaccard_Index_1step_in_bilateral, ...
    Male_Jaccard_Index_2step_in_out_bilateral, Male_Jaccard_Index_2step_in_bilateral];

herm_mf_data = [log(Herm_MFdiff_3node_struc_bilateral),log(Herm_MFdiff_3node_func_bilateral)...
    log(Herm_MFdiff_4node_struc_bilateral),log(Herm_MFdiff_4node_func_bilateral)];
male_mf_data = [log(Male_MFdiff_3node_struc_bilateral),log(Male_MFdiff_3node_func_bilateral)...
    log(Male_MFdiff_4node_struc_bilateral),log(Male_MFdiff_4node_func_bilateral)];

for i_redun=1:6
    if i_redun == 1
        data = herm_pr_data;
    elseif i_redun ==2
        data = male_pr_data;
    elseif i_redun ==3
        data = herm_cs_data;        
    elseif i_redun ==4
        data = male_cs_data;
    elseif i_redun ==5
        data = herm_mf_data;     
    elseif i_redun ==6
        data = male_mf_data;      
    end
X = normalize(data);

% PCA
[coeff, score, latent, tsquared, explained, mu] = pca(X);

target_sex = 'herm';
if strcmp(target_sex,'herm')
    idx_sensory = strcmp(Herm_celltype,'sensory neuron');
    idx_motor =  strcmp(Herm_celltype,'motorneuron');
    idx_inter =  strcmp(Herm_celltype,'interneuron');
    idx_sensory = idx_sensory(Herm_LRU(:,1)==1);
    idx_motor = idx_motor(Herm_LRU(:,1)==1);
    idx_inter = idx_inter(Herm_LRU(:,1)==1);
else
    idx_sensory = strcmp(Male_celltype,'sensory neuron');
    idx_motor =  strcmp(Male_celltype,'motorneuron');
    idx_inter =  strcmp(Male_celltype,'interneuron');
    idx_sensory = idx_sensory(Male_LRU(:,1)==1);
    idx_motor = idx_motor(Male_LRU(:,1)==1);
    idx_inter = idx_inter(Male_LRU(:,1)==1);
end


figure(1);
subplot(6,1,i_redun)
x=[1,2,3];
data=[mean(score(idx_sensory,1)),mean(score(idx_inter,1)),mean(score(idx_motor,1))];
b=bar(x,data);
b.FaceColor=[0.8,0.8,0.8];
hold on
err=[std(score(idx_sensory,1))/sqrt(length(idx_sensory)),std(score(idx_inter,1))/sqrt(length(idx_inter)),std(score(idx_motor,1))/sqrt(length(idx_motor))];
errlow=data-err;
errhigh=data+err;
er = errorbar(x,data,-err,err);
er.Color = [0 0 0];
er.LineStyle = 'none';
ylim([min(errlow)-0.1,max(errhigh)+0.1])


figure(2); subplot(3,2,i_redun);
 plotSpread(score(:,1),'categoryIdx',kmeans(score(:,1),3),'categoryColors',{'r','g','b'})
end

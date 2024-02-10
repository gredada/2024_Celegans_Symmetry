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


%% connectivity similarity (jaccard-index)
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

p_value_herm=custom_shuffle_test(Herm_Jaccard_Index_1step_in_all,Herm_Jaccard_Index_1step_in_bilateral,100000,1);
p_value_male=custom_shuffle_test(Male_Jaccard_Index_1step_in_all,Male_Jaccard_Index_1step_in_bilateral,100000,1);
disp(p_value_herm)
disp(p_value_male)

p_value_herm=custom_shuffle_test(Herm_Jaccard_Index_1step_in_out_all,Herm_Jaccard_Index_1step_in_out_bilateral,100000,1);
p_value_male=custom_shuffle_test(Male_Jaccard_Index_1step_in_out_all,Male_Jaccard_Index_1step_in_out_bilateral,100000,1);
disp(p_value_herm)
disp(p_value_male)

p_value_herm=custom_shuffle_test(Herm_Jaccard_Index_2step_in_all,Herm_Jaccard_Index_2step_in_bilateral,100000,1);
p_value_male=custom_shuffle_test(Male_Jaccard_Index_2step_in_all,Male_Jaccard_Index_2step_in_bilateral,100000,1);
disp(p_value_herm)
disp(p_value_male)

p_value_herm=custom_shuffle_test(Herm_Jaccard_Index_2step_in_out_all,Herm_Jaccard_Index_2step_in_out_bilateral,100000,1);
p_value_male=custom_shuffle_test(Male_Jaccard_Index_2step_in_out_all,Male_Jaccard_Index_2step_in_out_bilateral,100000,1);
disp(p_value_herm)
disp(p_value_male)
%% display for connectivity similarity

figure;
n=0.00001;
subplot(2,8,1); plotSpread({Herm_Jaccard_Index_1step_in_bilateral},'showMM',5,'distributionColor','k'); ylim([0,1])
subplot(2,8,2); plotSpread({Herm_Jaccard_Index_1step_in_all},'showMM',5,'distributionColor','k','binWidth',n); ylim([0,1])
subplot(2,8,3); plotSpread({Herm_Jaccard_Index_1step_in_out_bilateral},'showMM',5,'distributionColor','k'); ylim([0,1])
subplot(2,8,4); plotSpread({Herm_Jaccard_Index_1step_in_out_all},'showMM',5,'distributionColor','k','binWidth',n); ylim([0,1])

subplot(2,8,5); plotSpread({Herm_Jaccard_Index_2step_in_bilateral},'showMM',5,'distributionColor','k'); ylim([0,1])
subplot(2,8,6); plotSpread({Herm_Jaccard_Index_2step_in_all},'showMM',5,'distributionColor','k','binWidth',n); ylim([0,1])
subplot(2,8,7); plotSpread({Herm_Jaccard_Index_2step_in_out_bilateral},'showMM',5,'distributionColor','k'); ylim([0,1])
subplot(2,8,8); plotSpread({Herm_Jaccard_Index_2step_in_out_all},'showMM',5,'distributionColor','k','binWidth',n); ylim([0,1])


subplot(2,8,9); plotSpread({Male_Jaccard_Index_1step_in_bilateral},'showMM',5,'distributionColor','k'); ylim([0,1])
subplot(2,8,10); plotSpread({Male_Jaccard_Index_1step_in_all},'showMM',5,'distributionColor','k','binWidth',n); ylim([0,1])
subplot(2,8,11); plotSpread({Male_Jaccard_Index_1step_in_out_bilateral},'showMM',5,'distributionColor','k'); ylim([0,1])
subplot(2,8,12); plotSpread({Male_Jaccard_Index_1step_in_out_all},'showMM',5,'distributionColor','k','binWidth',n); ylim([0,1])

subplot(2,8,13); plotSpread({Male_Jaccard_Index_2step_in_bilateral},'showMM',5,'distributionColor','k'); ylim([0,1])
subplot(2,8,14); plotSpread({Male_Jaccard_Index_2step_in_all},'showMM',5,'distributionColor','k','binWidth',n); ylim([0,1])
subplot(2,8,15); plotSpread({Male_Jaccard_Index_2step_in_out_bilateral},'showMM',5,'distributionColor','k'); ylim([0,1])
subplot(2,8,16); plotSpread({Male_Jaccard_Index_2step_in_out_all},'showMM',5,'distributionColor','k','binWidth',n); ylim([0,1])

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

sum(mean(Herm_Jaccard_Index_1step_in_bilateral_null)>mean(Herm_Jaccard_Index_1step_in_bilateral))/1000
sum(mean(Herm_Jaccard_Index_1step_in_out_bilateral_null)>mean(Herm_Jaccard_Index_1step_in_out_bilateral))/1000
sum(mean(Herm_Jaccard_Index_2step_in_bilateral_null)>mean(Herm_Jaccard_Index_2step_in_bilateral))/1000
sum(mean(Herm_Jaccard_Index_2step_in_out_bilateral_null)>mean(Herm_Jaccard_Index_2step_in_out_bilateral))/1000
sum(mean(Male_Jaccard_Index_1step_in_bilateral_null)>mean(Male_Jaccard_Index_1step_in_bilateral))/1000
sum(mean(Male_Jaccard_Index_1step_in_out_bilateral_null)>mean(Male_Jaccard_Index_1step_in_out_bilateral))/1000
sum(mean(Male_Jaccard_Index_2step_in_bilateral_null)>mean(Male_Jaccard_Index_2step_in_bilateral))/1000
sum(mean(Male_Jaccard_Index_2step_in_out_bilateral_null)>mean(Male_Jaccard_Index_2step_in_out_bilateral))/1000
%% Symmetry/Asymmetry breaking for connectivity similarity
clear;clc;
cwd = pwd;
load('data/celegans_connectome.mat');

null_n=1000;
Herm_pair_n = sum(Herm_LRU(:,1));
Male_pair_n = sum(Male_LRU(:,1));

Herm_Jaccard_Index_1step_in_out_bilateral_asymlist  = zeros(Herm_pair_n,null_n,11);
Herm_Jaccard_Index_1step_in_bilateral_asymlist      = zeros(Herm_pair_n,null_n,11);
Male_Jaccard_Index_1step_in_out_bilateral_asymlist  = zeros(Male_pair_n,null_n,11);
Male_Jaccard_Index_1step_in_bilateral_asymlist      = zeros(Male_pair_n,null_n,11);
Herm_Jaccard_Index_2step_in_out_bilateral_asymlist  = zeros(Herm_pair_n,null_n,11);
Herm_Jaccard_Index_2step_in_bilateral_asymlist      = zeros(Herm_pair_n,null_n,11);
Male_Jaccard_Index_2step_in_out_bilateral_asymlist  = zeros(Male_pair_n,null_n,11);
Male_Jaccard_Index_2step_in_bilateral_asymlist      = zeros(Male_pair_n,null_n,11);

Herm_Jaccard_Index_1step_in_out_bilateral_symlist  = zeros(Herm_pair_n,null_n,11);
Herm_Jaccard_Index_1step_in_bilateral_symlist      = zeros(Herm_pair_n,null_n,11);
Male_Jaccard_Index_1step_in_out_bilateral_symlist  = zeros(Male_pair_n,null_n,11);
Male_Jaccard_Index_1step_in_bilateral_symlist      = zeros(Male_pair_n,null_n,11);
Herm_Jaccard_Index_2step_in_out_bilateral_symlist  = zeros(Herm_pair_n,null_n,11);
Herm_Jaccard_Index_2step_in_bilateral_symlist      = zeros(Herm_pair_n,null_n,11);
Male_Jaccard_Index_2step_in_out_bilateral_symlist  = zeros(Male_pair_n,null_n,11);
Male_Jaccard_Index_2step_in_bilateral_symlist      = zeros(Male_pair_n,null_n,11);
j=0;
tic;
for p_asym=0:0.1:1
    j=j+1;
    disp(j)
    toc;
    for i=1:null_n
        temp_A=randmio_dir_ratio(Herm_Adj, Herm_cont_ind, p_asym, 0, 1000);
        temp_B=randmio_dir_ratio(Male_Adj, Male_cont_ind, p_asym, 0, 1000);
        % jaccard index (1-step)
        [Herm_Jaccard_Index_1step_in_out_bilateral_asymlist(:,i,j), ~, ...
            Herm_Jaccard_Index_1step_in_bilateral_asymlist(:,i,j), ~] ...
            = asymtool_jaccard_1step(temp_A,Herm_LRU);
        [Male_Jaccard_Index_1step_in_out_bilateral_asymlist(:,i,j), ~, ...
            Male_Jaccard_Index_1step_in_bilateral_asymlist(:,i,j), ~] ...
            = asymtool_jaccard_1step(temp_B,Male_LRU);
        
        % jaccard index (2-step)
        [Herm_Jaccard_Index_2step_in_out_bilateral_asymlist(:,i,j), ~, ...
            Herm_Jaccard_Index_2step_in_bilateral_asymlist(:,i,j), ~] ...
            = asymtool_jaccard_2step(temp_A,Herm_LRU);
        [Male_Jaccard_Index_2step_in_out_bilateral_asymlist(:,i,j), ~, ...
            Male_Jaccard_Index_2step_in_bilateral_asymlist(:,i,j), ~] ...
            = asymtool_jaccard_2step(temp_B,Male_LRU);
    end
end
j=0;
tic;
for p_sym=0:0.1:1
    j=j+1;
    disp(j)
    toc;
    for i=1:null_n
        temp_A=randmio_dir_ratio(Herm_Adj, Herm_cont_ind, 0, p_sym, 1000);
        temp_B=randmio_dir_ratio(Male_Adj, Male_cont_ind, 0, p_sym, 1000);
        % jaccard index (1-step)
        [Herm_Jaccard_Index_1step_in_out_bilateral_symlist(:,i,j), ~, ...
            Herm_Jaccard_Index_1step_in_bilateral_symlist(:,i,j), ~] ...
            = asymtool_jaccard_1step(temp_A,Herm_LRU);
        [Male_Jaccard_Index_1step_in_out_bilateral_symlist(:,i,j), ~, ...
            Male_Jaccard_Index_1step_in_bilateral_symlist(:,i,j), ~] ...
            = asymtool_jaccard_1step(temp_B,Male_LRU);
        
        % jaccard index (2-step)
        [Herm_Jaccard_Index_2step_in_out_bilateral_symlist(:,i,j), ~, ...
            Herm_Jaccard_Index_2step_in_bilateral_symlist(:,i,j), ~] ...
            = asymtool_jaccard_2step(temp_A,Herm_LRU);
        [Male_Jaccard_Index_2step_in_out_bilateral_symlist(:,i,j), ~, ...
            Male_Jaccard_Index_2step_in_bilateral_symlist(:,i,j), ~] ...
            = asymtool_jaccard_2step(temp_B,Male_LRU);
    end
end

cd results
save jaccard_index_break *Jaccard_Index*
cd(cwd);
%%
figure;

data = squeeze(mean(Herm_Jaccard_Index_1step_in_bilateral_symlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,4,1);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')

data = squeeze(mean(Herm_Jaccard_Index_1step_in_out_bilateral_symlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,4,2);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')

data = squeeze(mean(Herm_Jaccard_Index_2step_in_bilateral_symlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,4,3);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')

data = squeeze(mean(Herm_Jaccard_Index_2step_in_out_bilateral_symlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,4,4);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')

data = squeeze(mean(Herm_Jaccard_Index_1step_in_bilateral_asymlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,4,5);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')

data = squeeze(mean(Herm_Jaccard_Index_1step_in_out_bilateral_asymlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,4,6);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')

data = squeeze(mean(Herm_Jaccard_Index_2step_in_bilateral_asymlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,4,7);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')

data = squeeze(mean(Herm_Jaccard_Index_2step_in_out_bilateral_asymlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,4,8);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')

%
figure;

data = squeeze(mean(Male_Jaccard_Index_1step_in_bilateral_symlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,4,1);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')

data = squeeze(mean(Male_Jaccard_Index_1step_in_out_bilateral_symlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,4,2);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')

data = squeeze(mean(Male_Jaccard_Index_2step_in_bilateral_symlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,4,3);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')

data = squeeze(mean(Male_Jaccard_Index_2step_in_out_bilateral_symlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,4,4);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')

data = squeeze(mean(Male_Jaccard_Index_1step_in_bilateral_asymlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,4,5);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')

data = squeeze(mean(Male_Jaccard_Index_1step_in_out_bilateral_asymlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,4,6);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')

data = squeeze(mean(Male_Jaccard_Index_2step_in_bilateral_asymlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,4,7);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')

data = squeeze(mean(Male_Jaccard_Index_2step_in_out_bilateral_asymlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,4,8);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')

%%
data = squeeze(mean(Male_Jaccard_Index_2step_in_bilateral_asymlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,2,3);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')
% xlim([-.1,1.1])
% ylim([.12,.19])
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
p_value_male=custom_shuffle_test(Male_MFdiff_3node_struc_all,Male_MFdiff_3node_struc_bilateral,100000,0);
disp(p_value_herm)
disp(p_value_male)

p_value_herm=custom_shuffle_test(Herm_MFdiff_3node_func_all,Herm_MFdiff_3node_func_bilateral,100000,0);
p_value_male=custom_shuffle_test(Male_MFdiff_3node_func_all,Male_MFdiff_3node_func_bilateral,100000,0);
disp(p_value_herm)
disp(p_value_male)

p_value_herm=custom_shuffle_test(Herm_MFdiff_4node_struc_all,Herm_MFdiff_4node_struc_bilateral,100000,0);
p_value_male=custom_shuffle_test(Male_MFdiff_4node_struc_all,Male_MFdiff_4node_struc_bilateral,100000,0);
disp(p_value_herm)
disp(p_value_male)

p_value_herm=custom_shuffle_test(Herm_MFdiff_4node_func_all,Herm_MFdiff_4node_func_bilateral,100000,0);
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

sum(mean(Herm_MFdiff_3node_struc_bilateral_null)>mean(Herm_MFdiff_3node_struc_bilateral))/1000
sum(mean(Herm_MFdiff_3node_func_bilateral_null)<mean(Herm_MFdiff_3node_func_bilateral))/1000
sum(mean(Herm_MFdiff_4node_struc_bilateral_null)>mean(Herm_MFdiff_4node_struc_bilateral))/1000
sum(mean(Herm_MFdiff_4node_func_bilateral_null)<mean(Herm_MFdiff_4node_func_bilateral))/1000
sum(mean(Male_MFdiff_3node_struc_bilateral_null)>mean(Male_MFdiff_3node_struc_bilateral))/1000
sum(mean(Male_MFdiff_3node_func_bilateral_null)<mean(Male_MFdiff_3node_func_bilateral))/1000
sum(mean(Male_MFdiff_4node_struc_bilateral_null)>mean(Male_MFdiff_4node_struc_bilateral))/1000
sum(mean(Male_MFdiff_4node_func_bilateral_null)<mean(Male_MFdiff_4node_func_bilateral))/1000

%% Symmetry/Asymmetry breaking for motif fingerprint difference
clear;clc;
cwd = pwd;
load('data/celegans_connectome.mat');

null_n=100;
Herm_pair_n = sum(Herm_LRU(:,1));
Male_pair_n = sum(Male_LRU(:,1));

% Herm_MFdiff_3node_struc_bilateral_asymlist = zeros(Herm_pair_n,null_n,11);
% Herm_MFdiff_3node_func_bilateral_asymlist = zeros(Herm_pair_n,null_n,11);
% Male_MFdiff_3node_struc_bilateral_asymlist = zeros(Male_pair_n,null_n,11);
% Male_MFdiff_3node_func_bilateral_asymlist = zeros(Male_pair_n,null_n,11);
Herm_MFdiff_4node_struc_bilateral_asymlist = zeros(Herm_pair_n,null_n,11);
Herm_MFdiff_4node_func_bilateral_asymlist = zeros(Herm_pair_n,null_n,11);
Male_MFdiff_4node_struc_bilateral_asymlist = zeros(Male_pair_n,null_n,11);
Male_MFdiff_4node_func_bilateral_asymlist = zeros(Male_pair_n,null_n,11);

% Herm_MFdiff_3node_struc_bilateral_symlist = zeros(Herm_pair_n,null_n,11);
% Herm_MFdiff_3node_func_bilateral_symlist = zeros(Herm_pair_n,null_n,11);
% Male_MFdiff_3node_struc_bilateral_symlist = zeros(Male_pair_n,null_n,11);
% Male_MFdiff_3node_func_bilateral_symlist = zeros(Male_pair_n,null_n,11);
Herm_MFdiff_4node_struc_bilateral_symlist = zeros(Herm_pair_n,null_n,11);
Herm_MFdiff_4node_func_bilateral_symlist = zeros(Herm_pair_n,null_n,11);
Male_MFdiff_4node_struc_bilateral_symlist = zeros(Male_pair_n,null_n,11);
Male_MFdiff_4node_func_bilateral_symlist = zeros(Male_pair_n,null_n,11);

tic;
j=0;
for p_asym=0:0.1:1
    j=j+1;
    disp(j)
    toc;
    
    for i=1:null_n
        
        null_Herm_Adj=randmio_dir_ratio(Herm_Adj, Herm_cont_ind, p_asym, 0, 1000);
        null_Male_Adj=randmio_dir_ratio(Male_Adj, Male_cont_ind, p_asym, 0, 1000);
        
        %         [Herm_MFdiff_3node_struc_bilateral_asymlist(:,i,j), ~, Herm_MFdiff_3node_func_bilateral_asymlist(:,i,j), ~] ...
        %             = asymtool_motif_fingerprint_difference(null_Herm_Adj,Herm_LRU,3,true,false,true,true);
        [Herm_MFdiff_4node_struc_bilateral_asymlist(:,i,j), ~, Herm_MFdiff_4node_func_bilateral_asymlist(:,i,j), ~] ...
            = asymtool_motif_fingerprint_difference(null_Herm_Adj,Herm_LRU,4,true,false,true,true);
        
        %         [Male_MFdiff_3node_struc_bilateral_asymlist(:,i,j), ~, Male_MFdiff_3node_func_bilateral_asymlist(:,i,j), ~] ...
        %             = asymtool_motif_fingerprint_difference(null_Male_Adj,Male_LRU,3,true,false,true,true);
        [Male_MFdiff_4node_struc_bilateral_asymlist(:,i,j), ~, Male_MFdiff_4node_func_bilateral_asymlist(:,i,j), ~] ...
            = asymtool_motif_fingerprint_difference(null_Male_Adj,Male_LRU,4,true,false,true,true);
        
    end
end
tic;
j=0;
for p_sym=0:0.1:1
    j=j+1;
    disp(j)
    toc;
    
    for i=1:null_n
        null_Herm_Adj=randmio_dir_ratio(Herm_Adj, Herm_cont_ind, 0, p_sym, 1000);
        null_Male_Adj=randmio_dir_ratio(Male_Adj, Male_cont_ind, 0, p_sym, 1000);
        
        %         [Herm_MFdiff_3node_struc_bilateral_symlist(:,i,j), ~, Herm_MFdiff_3node_func_bilateral_symlist(:,i,j), ~] ...
        %             = asymtool_motif_fingerprint_difference(null_Herm_Adj,Herm_LRU,3,true,false,true,true);
        [Herm_MFdiff_4node_struc_bilateral_symlist(:,i,j), ~, Herm_MFdiff_4node_func_bilateral_symlist(:,i,j), ~] ...
            = asymtool_motif_fingerprint_difference(null_Herm_Adj,Herm_LRU,4,true,false,true,true);
        
        %         [Male_MFdiff_3node_struc_bilateral_symlist(:,i,j), ~, Male_MFdiff_3node_func_bilateral_symlist(:,i,j), ~] ...
        %             = asymtool_motif_fingerprint_difference(null_Male_Adj,Male_LRU,3,true,false,true,true);
        [Male_MFdiff_4node_struc_bilateral_symlist(:,i,j), ~, Male_MFdiff_4node_func_bilateral_symlist(:,i,j), ~] ...
            = asymtool_motif_fingerprint_difference(null_Male_Adj,Male_LRU,4,true,false,true,true);
    end
end

cd results
save motif_fingerprint_difference_break *MFdiff* '-append'
cd(cwd);

data = squeeze(mean(Herm_MFdiff_3node_struc_bilateral_symlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,2,1);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')

data = squeeze(mean(Herm_MFdiff_3node_func_bilateral_symlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,2,2);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')

data = squeeze(mean(Herm_MFdiff_3node_struc_bilateral_asymlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,2,3);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')

data = squeeze(mean(Herm_MFdiff_3node_func_bilateral_asymlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,2,4);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')


data = squeeze(mean(Male_MFdiff_3node_struc_bilateral_symlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,2,1);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')

data = squeeze(mean(Male_MFdiff_3node_func_bilateral_symlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,2,2);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')

data = squeeze(mean(Male_MFdiff_3node_struc_bilateral_asymlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,2,3);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')

data = squeeze(mean(Male_MFdiff_3node_func_bilateral_asymlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,2,4);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')
%% pair-robustness
clear;clc;
cwd = pwd;
load('data/celegans_connectome.mat');

[Herm_Pair_Redundancy_bilateral,Herm_Pair_Redundancy_all] = asymtool_pair_redundancy_index2(Herm_Adj,Herm_LRU);
[Male_Pair_Redundancy_bilateral,Male_Pair_Redundancy_all] = asymtool_pair_redundancy_index2(Male_Adj,Male_LRU);

cd results
save pair_redundancy_index2 *Pair_Redundancy*
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
        = asymtool_pair_redundancy_index2(null_Herm_Adj,Herm_LRU,true,false);
    [Male_Pair_Redundancy_bilateral_null(:,i),~] ...
        = asymtool_pair_redundancy_index2(null_Male_Adj,Male_LRU,true,false);
end

cd results
save pair_redundancy_index_null2 *Pair_Redundancy*
cd(cwd);

%% Symmetry/Asymmetry breaking for pair-robustness
clear;clc;
cwd = pwd;
load('data/celegans_connectome.mat');

null_n=100;
Herm_pair_n = sum(Herm_LRU(:,1));
Male_pair_n = sum(Male_LRU(:,1));

Herm_Pair_Redundancy_bilateral_asymlist = zeros(Herm_pair_n,null_n,11);
Male_Pair_Redundancy_bilateral_asymlist = zeros(Male_pair_n,null_n,11);

Herm_Pair_Redundancy_bilateral_symlist = zeros(Herm_pair_n,null_n,11);
Male_Pair_Redundancy_bilateral_symlist = zeros(Male_pair_n,null_n,11);

tic;
j=0;
for p_asym=0:0.1:1
    j=j+1;
    disp(j)
    toc;
    for i=1:null_n
        null_Herm_Adj=randmio_dir_ratio(Herm_Adj, Herm_cont_ind, p_asym, 0, 1000);
        null_Male_Adj=randmio_dir_ratio(Male_Adj, Male_cont_ind, p_asym, 0, 1000);
        
        [Herm_Pair_Redundancy_bilateral_asymlist(:,i,j),~] ...
            = asymtool_pair_redundancy_index2(null_Herm_Adj,Herm_LRU,true,false);
        [Male_Pair_Redundancy_bilateral_asymlist(:,i,j),~] ...
            = asymtool_pair_redundancy_index2(null_Male_Adj,Male_LRU,true,false);
    end
end
j=0;
for p_sym=0:0.1:1
    j=j+1;
    disp(j)
    toc;
    for i=1:null_n
        null_Herm_Adj=randmio_dir_ratio(Herm_Adj, Herm_cont_ind, 0, p_sym, 1000);
        null_Male_Adj=randmio_dir_ratio(Male_Adj, Male_cont_ind, 0, p_sym, 1000);
        
        [Herm_Pair_Redundancy_bilateral_symlist(:,i,j),~] ...
            = asymtool_pair_redundancy_index2(null_Herm_Adj,Herm_LRU,true,false);
        [Male_Pair_Redundancy_bilateral_symlist(:,i,j),~] ...
            = asymtool_pair_redundancy_index2(null_Male_Adj,Male_LRU,true,false);
    end
end

cd results
save pair_redundancy_index_break *Pair_Redundancy*
cd(cwd);

figure;

data = squeeze(mean(Herm_Pair_Redundancy_bilateral_symlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,2,1);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')

data = squeeze(mean(Male_Pair_Redundancy_bilateral_symlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,2,2);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')

data = squeeze(mean(Herm_Pair_Redundancy_bilateral_asymlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,2,3);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')

data = squeeze(mean(Male_Pair_Redundancy_bilateral_asymlist));
errmin=mink(data,1)-mean(data);
errmax=maxk(data,1)-mean(data);
subplot(2,2,4);errorbar((0:0.1:1),mean(data),errmin,errmax,'-k.')

%% display for pair-robustness

figure; plotSpread({Herm_Pair_Redundancy_bilateral,Herm_Pair_Redundancy_all,[],Male_Pair_Redundancy_bilateral,Male_Pair_Redundancy_all},'showMM',5,'distributionColor','k')
sum(mean(Herm_Pair_Redundancy_bilateral_null)>mean(Herm_Pair_Redundancy_bilateral))/1000
sum(mean(Male_Pair_Redundancy_bilateral_null)>mean(Male_Pair_Redundancy_bilateral))/1000

label=Herm_label;
LRU=Herm_LRU;
temp=label(find(LRU(:,1)));
for i=1:length(temp)
    temp{i}=temp{i}(1:end-1);
end
find(strcmp(temp,'AVA'))

figure;
temp_herm=idx_herm;
% temp_herm=ones(length(idx_herm),1);
temp_herm(14)=4;temp_herm(31)=5;
temp_herm(21)=1;temp_herm(22)=1;
temp_male=idx_male;
% temp_male=ones(length(idx_male),1);
temp_male(14)=4;temp_male(31)=4;
temp_male(21)=1;temp_male(22)=5;
plotSpread([{score_herm} {score_male}],'categoryIdx',[temp_herm;temp_male],'categoryColors',{'r','g','b','b','b'},'categoryMarkers',{'o','.','.','^','^'},'xNames',{'Herm','Male'});
%% PCA group analysis
clear; clc;
close all

% data
load('data/celegans_connectome.mat');

load('results/pair_redundancy_index2.mat');
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
    if i_redun == 5
        data = herm_pr_data;
        target_sex = 'herm';
    elseif i_redun ==6
        data = male_pr_data;
        target_sex = 'male';
        
    elseif i_redun ==1
        data = herm_cs_data;
        target_sex = 'herm';
        
    elseif i_redun ==2
        data = male_cs_data;
        target_sex = 'male';
        
    elseif i_redun ==3
        data = herm_mf_data;
        target_sex = 'herm';
        
    elseif i_redun ==4
        data = male_mf_data;
        target_sex = 'male';
        
    end
    X = normalize(data);
    
    % PCA
    [coeff, score, latent, tsquared, explained, mu] = pca(X);
    
    if strcmp(target_sex,'herm')
        idx_sensory = strcmp(Herm_celltype,'sensory neuron');
        idx_motor =  strcmp(Herm_celltype,'motorneuron');
        idx_inter =  strcmp(Herm_celltype,'interneuron');
        idx_sensory = idx_sensory(Herm_LRU(:,1)==1);
        idx_motor = idx_motor(Herm_LRU(:,1)==1);
        idx_inter = idx_inter(Herm_LRU(:,1)==1);
        label=Herm_label;
        LRU=Herm_LRU;
    else
        idx_sensory = strcmp(Male_celltype,'sensory neuron');
        idx_motor =  strcmp(Male_celltype,'motorneuron');
        idx_inter =  strcmp(Male_celltype,'interneuron');
        idx_sensory = idx_sensory(Male_LRU(:,1)==1);
        idx_motor = idx_motor(Male_LRU(:,1)==1);
        idx_inter = idx_inter(Male_LRU(:,1)==1);
        label=Male_label;
        LRU=Male_LRU;
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
    
    
    %
    % alllow=1;
    % x_s=score(idx_sensory,1);
    % x_i=score(idx_inter,1);
    % x_m=score(idx_motor,1);
    % [~,~,~,stats] = ttest2(x_s,x_i);
    % t_obs1 = stats.tstat;
    % [~,~,~,stats] = ttest2(x_s,x_m);
    % t_obs2 = stats.tstat;
    % [~,~,~,stats] = ttest2(x_m,x_i);
    % t_obs3 = stats.tstat;
    % t_obs=max([abs(t_obs1),abs(t_obs2),abs(t_obs3)]);
    %
    % n=100000;
    % parfor i=1:n
    %     perm_ind=randperm(length(idx_sensory));
    %
    %     y_s=score(idx_sensory(perm_ind),1);
    %     y_i=score(idx_inter(perm_ind),1);
    %     y_m=score(idx_motor(perm_ind),1);
    %
    %
    %     [~,~,~,stats] = ttest2(y_s,y_i);
    %     t_tar1 = stats.tstat;
    %     [~,~,~,stats] = ttest2(y_s,y_m);
    %     t_tar2 = stats.tstat;
    %     [~,~,~,stats] = ttest2(y_m,y_i);
    %     t_tar3 = stats.tstat;
    %     t_tar=max([abs(t_tar1),abs(t_tar2),abs(t_tar3)]);
    %
    %     temp_t(i)= abs(t_tar);
    %
    %
    % end
    %
    % p_val_1=sum(temp_t >= abs(t_obs1))/n
    % p_val_2=sum(temp_t >= abs(t_obs2))/n
    % p_val_3=sum(temp_t >= abs(t_obs3))/n
    %
    %
    % temp_p(i_redun,:)=[p_val_1,p_val_2,p_val_3];
    rng(123);
    temp_kmeans=kmeans(score(:,1),3);
    [~,temp_red_idx]=sort([mean(score(temp_kmeans==1,1)),mean(score(temp_kmeans==2,1)),mean(score(temp_kmeans==3,1))]);
    temp=label(find(LRU(:,1)));
    for i=1:length(temp)
        temp{i}=temp{i}(1:end-1);
    end
    group{1,i_redun}=temp(temp_kmeans==temp_red_idx(1));
    group{2,i_redun}=temp(temp_kmeans==temp_red_idx(2));
    group{3,i_redun}=temp(temp_kmeans==temp_red_idx(3));
    
    
    sorted_idx = zeros(size(temp_kmeans));
    for i = 1:3
        sorted_idx(temp_kmeans == temp_red_idx(i)) = i;
    end
    %         sorted_idx(14)=4;
    %         sorted_idx(31)=5;
    %         sorted_idx(21)=6;
    %         sorted_idx(22)=7;
    %         if i_redun==5
    %             sorted_idx(22)=3;
    %
    %         end
    figure(2); subplot(3,2,i_redun);
    if i_redun == 3
        plotSpread(score(:,1),'categoryIdx',sorted_idx,'categoryColors',{'r','g','b'})
                    ylim([-5,5])

    elseif i_redun == 4
        plotSpread(score(:,1),'categoryIdx',sorted_idx,'categoryColors',{'r','g','b'})
                    ylim([-5,5])

    elseif i_redun == 5
        plotSpread(score(:,1),'categoryIdx',sorted_idx,'categoryColors',{'b','g','r'})
                    ylim([-2,6])

    elseif i_redun == 1
        plotSpread(score(:,1),'categoryIdx',sorted_idx,'categoryColors',{'b','g','r'})
            ylim([-5,5])

    elseif i_redun == 2
        plotSpread(score(:,1),'categoryIdx',sorted_idx,'categoryColors',{'b','g','r'})
                    ylim([-5,5])

    elseif i_redun == 6
        plotSpread(score(:,1),'categoryIdx',sorted_idx,'categoryColors',{'b','g','r'})
                            ylim([-2,6])

    else
        
    end
    %  [~,~,temp]=anova1(score,idx_inter+2*idx_sensory+3*idx_motor)
    % c=multcompare(temp)
    [score(14),score(31),score(21),score(22)]
    % temp_herm(14)=4;temp_herm(31)=5;
    % temp_herm(21)=1;temp_herm(22)=1;
    % temp_male=idx_male;
    % % temp_male=ones(length(idx_male),1);
    % temp_male(14)=4;temp_male(31)=4;
    % temp_male(21)=1;temp_male(22)=5;
    
    
    
end

% group
clear 
close all;

addpath(genpath('/nobackup/dcarrol2/MATLAB'));

modelDir = '/nobackup/dcarrol2/v05_latest/darwin3/run/diags/monthly/';

codeDir = '/nobackup/dcarrol2/evaluation/m_files/model_obs_compare/BGC-Argo/';

%% 

cd(codeDir);

model_obs_compare_profiles_NAS(1,modelDir);

%%

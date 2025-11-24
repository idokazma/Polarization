addpath(genpath(pwd))
clear all


params.multisources = false;

params.only_center = true;
single_scatterer_main;
disp(1)

clear all
params.multisources = false;

params.only_center = false;
single_scatterer_main;
disp(2)





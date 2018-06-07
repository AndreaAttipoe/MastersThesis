clear all; clc; close all;
%% DimensionalityAnalysis - Andrea ATTIPOE - Master's Thesis 2017-2018.
% Runs the MSD_2D3D function Nwalk times and computes the diffusion 
% coefficients averages and stds.

%% Parameters
Nwalks=10^2;
Nsteps=10^3;
Navg=10;
D_2Ds=zeros(Nwalks,1);
D_3Ds=zeros(Nwalks,1);
D_2Deulers=zeros(Nwalks,1);
%% MSDs computation
for i=1:Nwalks
    i
    [D_2Ds(i), D_3Ds(i), D_2Deulers(i)]=MSD2D3D(Nsteps,Navg,0);
end
%% Means
mean2D=mean(D_2Ds)
mean3D=mean(D_3Ds)
mean2Deuler=mean(D_2Deulers)
std2D=std(D_2Ds)
std3D=std(D_3Ds)
std2Deuler=std(D_2Deulers)

save('normDimMean1StdHalf.mat');
clc; clear all; close all;
%% MicrorheologyMSD - Andrea ATTIPOE - Master's Thesis 2017-2018.
% Imports the files containing the MSD data (they should be in the 'Data' folder
% at the root of the script), plots the linear and loglog MSD curves, computes
% the logarithmic slope coefficient, the Diffusion coefficient and the dynamic
% viscosity coefficient with its associated error.
% The files provided in Data are related to the calibration with glycerol.
cd('Data');

% Parameters
nFPS = 100; %[fps]
umPixel = 0.0446; % [micro m/pixels] Calibration constant
T=273.15+20; %[K]
limitMSD=100; %[frames] Length in frames over which the regression is performed

%% Importing Recording
% File reading
[filename, path_n] = uigetfile('*.txt', 'MultiSelect', 'on');
if ~iscell(filename)
    filename1 = {filename};
end
for j = 1:size(filename,2);
    data = dlmread(filename{j},'\t',2,0);
    Time{j} = data(:,1)/nFPS; % Scaling
    MSD{j} = data(:,3:end)*(umPixel^2); % Scaling
end

%% MSD Plots
% Logarithmic plot
figure1 = figure;
axes1= axes('Parent', figure1);
set(gcf,'Units','centimeters');
set(gcf,'Position',[0.0 0.0 60 60*3/4]);
set(gcf,'PaperPosition',[0.0 0.0 60 60*3/4]);
grid on;
box on;
set(gca,'Fontsize',24);
% Colors
colors=hsv(size(MSD,2));
hold on;
loglog((1/nFPS):0.1:Time{j}(end), (1/nFPS):0.1:Time{j}(end), 'k', ...
'LineWidth',1.5);
for j = 1:size(filename,2);
    loglog(Time{j},MSD{j},'o', 'color', colors(j,:), 'MarkerSize', 3);
end
hold off;
xlabel('Time $\tau$ [s]','Interpreter','latex');
ylabel('MSD [$\mu$m$^2$]','Interpreter','latex');
lgd=legend('Newtonian', 'Data','Location','best');
set(lgd,'FontSize',15);
title('Logarithmic plot', 'Interpreter', 'latex');
% Set the remaining axes properties
set(axes1,'XGrid','on','XMinorTick','on','XScale','log','YGrid','on', ...
'YMinorTick','on','YScale','log');

% Linear plot
figure2 = figure;
axes2= axes('Parent', figure2);
set(gcf,'Units','centimeters');
set(gcf,'Position',[0.0 0.0 60 60*3/4]);
set(gcf,'PaperPosition',[0.0 0.0 60 60*3/4]);
grid on;
box on;
set(gca,'Fontsize',24);
% Colors
colors=hsv(size(MSD,2));
hold on;
for j = 1:size(filename,2);
    plot(Time{j},MSD{j},'o', 'color', colors(j,:), 'MarkerSize', 3);
end
hold off;
xlabel('Time $\tau$ [s]','Interpreter','latex');
ylabel('MSD [$\mu$m$^2$]','Interpreter','latex');
set(lgd,'FontSize',15);
title('Linear plot', 'Interpreter', 'latex');
% Set the remaining axes properties
set(axes2,'XGrid','on','XMinorTick','on','XScale','lin','YGrid','on', ...
'YMinorTick','on','YScale','lin');

%% Average and slope computation
meanMSD=zeros(limitMSD,1);
for j=1:size(filename,2);
    meanMSD=meanMSD+MSD{j}(1:limitMSD,1);
end
meanMSD=meanMSD/size(filename,2);
[MeanFit,data]=polyfit(Time{j}(1:limitMSD),(10^-12)*meanMSD,1);
fprintf('Coefficients : Slope = %f, Offset = %f, Norm of residuals : %f\n', ...
MeanFit(1),MeanFit(2),data.normr);
meanSlope=MeanFit(1);

%% Log slope
[LogFit,dataLog]=polyfit(log(Time{1}(1:limitMSD)),log(meanMSD),1);
fprintf('Logarithmic slope = %f, Norm of residuals : %f\n',LogFit(1), ...
dataLog.normr);

%% Dynamic viscosity coefficient
kB=1.38064852*10^(-23); % [m^2 kg/s^2.K] Boltzman constant
R=(0.74/2)*10^(-6); %[m]
D=(meanSlope)/4 %[m^2/s]
eta=kB*T/(6*pi*R*D)*1000 %[mPa.s]

%% Error computation
deltaR=0.0025*10^-6; %[m]
deltaT=0.5; %[K]
regressedMSD=MeanFit(1)*Time{1}(1:limitMSD)+MeanFit(2);
    for j = 1:size(filename,2);
        errorMSD((1+(j-1)*limitMSD):(limitMSD+(j-1)*limitMSD))=10^(-12) ...
        *MSD{j}(1:limitMSD)-regressedMSD;
    end
meanErrorMSD=mean(errorMSD);
deltaD=std(errorMSD)/(4*sqrt(length(errorMSD))); %[m^2/s]
deltaEta=(abs(-kB*T/(6*pi*R*(D^2)))*deltaD+abs(-kB*T/(6*pi*(R^2)*D))*deltaR ...
+abs(kB/(6*pi*R*D))*deltaT)*1000 %[mPa.s]

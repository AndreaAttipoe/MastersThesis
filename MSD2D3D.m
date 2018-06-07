function [D_2D, D_3D, D_2Deuler] = MSD2D3D(Njumps, MSDlimit, plotter)
%% MSD2D3D - Andrea ATTIPOE - Master's Thesis 2017-2018.
% Computes the MSD of a random walk in a projected 2D, 3D and Euler angles
% projected 2D context. Njumps is the number of steps in the walk,
% MSDlimit the number of sampling times considered in the average,
% plotter determines if the random walk is displayed on a 3D plot.

%% Parameters
dr=1; %[l.u.]
dt=1; %[t.u.]
X=zeros(Njumps,1);
Y=zeros(Njumps,1);
Z=zeros(Njumps,1);
% Euler angles:
a=degtorad(5); %[rad]
b=degtorad(5); %[rad]
%% Random walk

theta=(2*pi)*rand(Njumps-1,1); %[rad]
phi=(pi)*rand(Njumps-1,1); %[rad]

if(plotter==1)
    figure;
end

for i=2:Njumps
    % Jumps
    X(i)=X(i-1)+dr*cos(theta(i-1))*sin(phi(i-1));
    Y(i)=Y(i-1)+dr*sin(theta(i-1))*sin(phi(i-1));
    Z(i)=Z(i-1)+dr*cos(phi(i-1));
    % Euler Plane Projection
    Xeuler(i)=X(i)*cos(b)+Z(i)*sin(b);
    Yeuler(i)=-X(i)*sin(a)*sin(b)+Y(i)*cos(a)+Z(i)*sin(a)*cos(b);
    Zeuler(i)=-X(i)*cos(a)*sin(b)-Y(i)*sin(a)+Z(i)*cos(a)*cos(b);
    if(plotter==1)
        % 3D Plot
        plot3(X(i),Y(i),Z(i),'.','Markersize',10,'MarkerEdgeColor','r')
        hold on
        line([X(i-1) X(i)], [Y(i-1) Y(i)],[Z(i-1) Z(i)],'linewidth',1);
        axis([ -10 10 -10 10 -10 10])
        view([-15 15])
        grid on
        xlabel('X [l.u.]','Interpreter','latex');
        ylabel('Y [l.u.]','Interpreter','latex');
        zlabel('Z [l.u.]','Interpreter','latex');
        set(gca,'Fontsize',18);
        Frame=getframe;
    end
end
if(plotter==1)
    movie(Frame);
end

%% MSD 2D
Tau=0;
for i=1:(Njumps-1)
    Tau=Tau+dt;
    time2D(i)=Tau;
    k=1;
    for j=1:Njumps
        if(j+Tau/dt<length(X))
            SD_2D(k)=(X(j+Tau/dt)-X(j))^2+(Y(j+Tau/dt)-Y(j))^2;
            k=k+1;
        end
    MSD_2D(i)=mean(SD_2D);
    end
end
%% MSD 3D
Tau=0;
for i=1:(Njumps-1)
    Tau=Tau+dt;
    time3D(i)=Tau;
    k=1;
    for j=1:Njumps
        if(j+Tau/dt<length(X))
            SD_3D(k)=(X(j+Tau/dt)-X(j))^2+(Y(j+Tau/dt)-Y(j))^2+(Z(j+Tau/dt) ...
            -Z(j))^2;
            k=k+1;
        end
    MSD_3D(i)=mean(SD_3D);
    end
end
%% MSD - Euler Plane
Tau=0;
for i=1:(Njumps-1)
    Tau=Tau+dt;
    time2Deuler(i)=Tau;
    k=1;
    for j=1:Njumps
        if(j+Tau/dt<length(Xeuler))
            SD_2Deuler(k)=(Xeuler(j+Tau/dt)-Xeuler(j))^2+(Yeuler(j+Tau/dt) ...
            -Yeuler(j))^2;
            k=k+1;
        end
    MSD_2Deuler(i)=mean(SD_2Deuler);
    end
end
%% D computation
[LinFit_2D,data]=polyfit(time2D(1:MSDlimit),MSD_2D(1:MSDlimit),1);
fprintf('2D Coefficients : a = %f, b = %f, Norm of residuals : %f\n', ...
LinFit_2D(1),LinFit_2D(2),data.normr);
[LinFit_3D,data]=polyfit(time3D(1:MSDlimit),MSD_3D(1:MSDlimit),1);
fprintf('3D Coefficients : a = %f, b = %f, Norm of residuals : %f\n', ...
LinFit_3D(1),LinFit_3D(2),data.normr);
[LinFit_2Deuler,data]=polyfit(time2Deuler(1:MSDlimit),MSD_2Deuler(1:MSDlimit),1);
fprintf('2D Euler Coefficients : a = %f, b = %f, Norm of residuals : %f\n', ...
LinFit_2Deuler(1),LinFit_2Deuler(2),data.normr);
D_2D=LinFit_2D(1)/4
D_3D=LinFit_3D(1)/6
D_2Deuler=LinFit_2Deuler(1)/4

if(plotter==1)
%% Log plot
figure1 = figure;
axes1= axes('Parent', figure1);
set(gcf,'Units','centimeters');
set(gcf,'Position',[0.0 0.0 60 60*3/4]);
set(gcf,'PaperPosition',[0.0 0.0 60 60*3/4]);
grid on;
box on;
set(gca,'Fontsize',34);
% Colors
hold on;
loglog(time2D,MSD_2D,'ro', 'MarkerSize', 4);
loglog(time3D,MSD_3D,'bo', 'MarkerSize', 4);
loglog(time2Deuler,MSD_2Deuler,'mo', 'MarkerSize', 4);
hold off;
xlabel('Time [t.u.]','Interpreter','latex');
ylabel('MSD [l.u.$^2$]','Interpreter','latex');
lgd=legend('2D', '3D', '2D Euler Projection', 'Location','best');
set(lgd,'FontSize',14);
title('Logarithmic plot', 'Interpreter', 'latex');
% Set the remaining axes properties
set(axes1,'XGrid','on','XMinorTick','on','XScale','log','YGrid','on',...
'YMinorTick','on','YScale','log');
end

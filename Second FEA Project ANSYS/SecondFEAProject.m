clear all
clc
close all

x = (0:0.1:1.5)

%GIVEN 

q = 1000; %N/m
L = 1.5; % m
r = 0.1; %m
I = (pi*r^4)/4; %m^4
E = 207*10^9; %Pa


v = (-q/(24*E*I))*(x.^4-4*L*x.^3+6*L^2*x.^2);
theta = (-q/(6*E*I))*(x.^3-3*L*x.^2+3*L^2*x);
Moment = (-q/2)*(x.^2-2*L*x+L^2);
Shear =  q*(x-L);

set(groot, 'DefaultTextInterpreter', 'LaTeX', ...
           'DefaultAxesTickLabelInterpreter', 'LaTeX', ...
           'DefaultAxesFontName', 'LaTeX', ...
           'DefaultLegendInterpreter', 'LaTeX', ...
           'defaultFigureColor','w');
       
fig1=figure; hold on; grid on; set(gca,'FontSize',14);


subplot(4,1,1);
plot(x,v,'LineWidth',3)
title('Displacement')
xlabel('x-dir (m)')
ylabel('y-dir (m)')

subplot(4,1,2);
plot(x,theta,'LineWidth',3)
title('Theta')
xlabel('x-dir (m)')
ylabel('Theta (rad)')

subplot(4,1,3);
plot(x,Moment,'LineWidth',3)
title('Moment')
xlabel('x-dir (m)')
ylabel('Moment (Nm)')

subplot(4,1,4);
plot(x,Shear,'LineWidth',3)
title('Shear')
xlabel('x-dir (m)')
ylabel('Shear (N)')

set(groot, 'Default', struct()) 
%% QUIZ 19 -- LUIS KLIGMAN

clear all
close all
clc

%% Create time interval & data table
t = linspace(eps,0.1,100);
x = linspace(eps,1,100);
[xtMesh,txMesh]=meshgrid(x,t);

%% New function of temperature vs time
y1 = @(x,t) 90 + exp(-pi^2*t).*cos(pi*x) + exp(-4*pi^2*t).*cos(2*pi*x);

%% Plot functiion
figure(1)
surf(xtMesh,txMesh, y1(xtMesh,txMesh),'FaceColor','interp',...
'EdgeColor','none',...
'FaceLighting','gouraud')
hold on
grid on
xlabel('x')
ylabel('time')
zlabel('temperature')
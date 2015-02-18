% Based on Surdutovich, Solov'yov (2010)
% Last edited by Ping Lin/Gabriel, 13th Feb
% Command line file to solve and plot pressure against wavefront radius

close all; clear all; clc; clear mem;

beta=0.86;
epsilon=0.9; % keV/nm
rho1=1; % g cm^-3
y=1.222; % gamma

t=0:2e-12:2e-10; % s
Radius=R(beta,epsilon,rho1,t); % wavefront radius in mm
pressure=pressure2(y,t);

figure(1)
plot(Radius*1e6,pressure*1e-11); % fig2 of paper
xlabel('R (nm)');
ylabel('P');

t_1=1e-12;
r=0:0.03:3; % radial distance from axis in nm
xi=xi(r,t_1,beta,rho1,epsilon); 
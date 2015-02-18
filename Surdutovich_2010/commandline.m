% Based on Surdutovich, Solov'yov (2010)
% Last edited by Ping Lin/Gabriel, 13th Feb
% Command line file to solve and plot pressure against wavefront radius

close all; clear all; clc; clear mem;

beta=0.86;
epsilon=0.9*1e3*1.602e-19*1e9; % keV/nm -> J/m
rho1=1*1e-3*1e6; % g cm^-3 -> kg m^-3
y=1.222; % gamma

t=0:2e-13:2e-10; % s
Radius=R(beta,epsilon,rho1,t); % wavefront radius in m
pressure=pressure2(y,t);

figure(1)
plot(Radius*1e9,pressure*1e-9); % fig2 of paper
xlabel('R (nm)');
ylabel('P (nN nm^-2)');

t_1=1e-12;
r=0:0.03e-9:3e-9; % radial distance from axis in m
xi=xi(r,t_1,beta,rho1,epsilon); 
G_1=(y+1)/(y-1); % G at t=1ps
Z_1=(2*y*(y-1))/(y+1)^2; % Z at t=1ps
pressure_2=pressure_r(xi,G_1,Z_1,t_1); % at t=1ps

figure(2)
plot(r*1e9,pressure_2*1e-9);
xlabel('r (nm)');
ylabel('P');
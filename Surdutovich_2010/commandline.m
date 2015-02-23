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
plot(Radius*1e9,pressure); % fig2 of paper
xlabel('R (nm)');
ylabel('P (nN nm^-2)');

t_1=1e-12;
r=0:0.03e-9:3e-9; % radial distance from axis in m
xi=xi(r,t_1,beta,rho1,epsilon); 

% G_1=(y+1)/(y-1); % G at t=1ps
% Z_1=(2*y*(y-1))/(y+1)^2; % Z at t=1ps
% pressure_2=pressure_r(xi,G_1,Z_1,t_1); % at t=1ps

% figure(2)
% plot(r*1e9,pressure_2*1e-9);
% xlabel('r (nm)');
% ylabel('P');

% V=1/y:0.01:2/y;
% xi2=xi_square(y,V);
% figure(3)
% plot(sqrt(xi2),V); 

% V=@(xilist) 1/y+xilist.^(2*y/(y-1));
% V_list=V(xi);
% G=@(xilist) xilist.^(2/(y-1));
% G_list=G(xi);
% Z=@(vlist) (y.*(y-1).*(1-vlist).*vlist.^2)./(2*(vlist.*y-1));
% Z_list=Z(V_list);
% pressure_2=pressure_r(xi,G_list,Z_list,t_1);
% 
% figure(3);
% plot(r*1e9,pressure_2*1e-9);
% xlabel('r (nm)');
% ylabel('P');

V_list2=get_V(y,xi);
G2=@(vlist) ((y-1)/(y+1)) * (((y-1)*(y.*vlist-1))./(y+1)).^(1/y) .* (((y-1).*(2-y.*vlist))./(2.*(1-vlist))).^(2/(2-y));
G_list2=G2(V_list2);
Z2=@(vlist) (y.*(y-1).*(1-vlist).*vlist.^2)./(2*(vlist.*y-1));
Z_list2=Z2(V_list2);
pressure_2=pressure_r(xi,G_list2,Z_list2,1);

figure(4);
plot(r*1e9,pressure_2*1e3); % why 4e3?
xlabel('r (nm)');
ylabel('P');

figure(5);
plot(r*1e9,gradient(pressure_2*1.2e5)); % why 1.2e5?
xlabel('r (nm)');
ylabel('dP/dr');
%-----------------------------------
% IT'S NOT WORKING!!!!!!!!!!11111
%-----------------------------------
% Based on Bortfeld 1997: An analytical approximation of the Bragg curve
% for therapeutic proton beams
% Last edited: 17th Jan
% Command line file to adjust parameters, plot graphs from functions

close all; clear all; clc; clear mem;

% declaring constants/parameters
E0=150; %E0 in units of MeV
alpha=2.2e-3;
p=1.77;
rho=1;
R0=range(alpha,E0,p); %R in units of cm
beta=0.012;
sigma=0.012*(R0^0.935);
gamma=0.6;
phi0=1000;
epsilon=0.2;

% testing functions
R=0:0.1:R0*1.1;
d=0:0.1:R0*1.1;
% R=0:0.1:R0;
% d=0:0.1:R0;
z=zetafunc(R0,d,sigma);

testmz=parafunc(-1/p,-z);
testmz2=parafunc2(-1/p,-z);
coeff = 2^(1/(2*p))*exp(-z.^2/4);
figure(1);
[ax,p1,p2]=plotyy(d,testmz,d,kummerU(-1/(2*p),1/2,0.5*z.^2),'plot','plot'); hold on;
ylabel(ax(1),'Parabolic function') % label left y-axis
ylabel(ax(2),'Kummer U Function') % label right y-axis
plot(d,coeff);
figure(2);
plot(d,testmz); hold on;
plot(d,testmz2);


 
D_z=dose(phi0,sigma,beta,alpha,gamma,E0,p,d,rho);
D_h20=doseh20(phi0,sigma,alpha,epsilon,E0,p,d);

figure(1);
plot(d,D_z);
ylabel('Dose per fluence');


figure(2);
plot(d,D_h20);

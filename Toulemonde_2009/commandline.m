% Based on Toulemonde, Surdutovich, Solov'yov (2009)
% Last edited by Ping Lin/Gabriel, 27th Jan
% Command line file to solve and plot temperature against time

close all; clear all; clc; clear mem;

% electronic subsystem
% t=0:0.001:0.1;
% initial_Te=310;
% [Te,t]=ode45(@electrontempfunc, t, initial_Te);
% figure(1);
% plot(t,Te);

% molecular subsystem
% t=0:0.001:0.1;
% initial_T=310;
% [T,t]=ode45(@moleculartempfunc, t, initial_T);
% figure(2);
% plot(t,T);

r=0:0.1:1000;
e1=energydensity(r,0.075,150,1);
figure(3)
semilogx(r,e1);
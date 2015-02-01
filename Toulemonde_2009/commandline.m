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

r=logspace(-7,0,250); % mm
[e1,e1_1]=energydensity(r,0.010,1,1);
[e10,e10_1]=energydensity(r,0.010,10,1);
[e20,e20_1]=energydensity(r,0.010,20,1);
[e50,e50_1]=energydensity(r,0.010,50,1);
[e100,e100_1]=energydensity(r,0.010,100,1);
figure(3);
loglog(r*1e6,e1,'b'); hold on; loglog(r*1e6,e1_1,'b--');
loglog(r*1e6,e10,'r'); loglog(r*1e6,e10_1,'r--');
loglog(r*1e6,e20,'Color',[0,0.5,1]); loglog(r*1e6,e20_1,'Color',[0,0.5,1],'LineStyle','--');
loglog(r*1e6,e50,'k'); loglog(r*1e6,e50_1,'k--');
loglog(r*1e6,e100,'m'); loglog(r*1e6,e100_1,'m--');
legend('1 MeV','1 MeV(without correction)','10 MeV','10 MeV(without correction)','20 MeV','20 MeV(without correction)','50 MeV','50 MeV(without correction)','100 MeV','100MeV(without correction)');
ylim([1e-8, 1e7]);

% [p,e,t] = initmesh('squareg'); 
% [p,e,t] = refinemesh('squareg',p,e,t); 
% u0 = zeros(2,1);  
% tlist = linspace(0,0.01,2); 
% bc=0; % boundary conditions
% c=[2;0;2;6e-3;0;6e-3]; % [Ke,K] in units of W K^-1 cm^-1 
% g=1; % electron-phonon coupling
% a=-1*g*ones(2,2);
% for i=1:2
%     a(i,i)=g;
% end
% d=zeros(2,2);
% d(1,1)=1; % Ce, in J K^-1 g^-1
% d(2,2)=4.18; % rho*C, in g cm^-3 J g^-1 K^-1
% r=0:0.1:1000;
% B=energydensity(r,0.075,150,1);
% f=[B;0];
% u1 = parabolic(u0,tlist,'bc',p,e,t,c,a,f,d);
% Based on Bortfeld 1997: An analytical approximation of the Bragg curve
% for therapeutic proton beams
% Last edited: 22nd Jan
% Command line file to adjust parameters, plot graphs from functions

close all; clear all; clc; clear mem;

% Declaring constants/parameters
E0=150; %E0 in units of MeV
alpha=2.2e-3;
p=1.77;
rho=1; %mass density of medium, g/cm^3
R0=range(alpha,E0,p); %R in units of cm
beta=0.012;
gamma=0.6; %fraction of energy released in the nonelastic nuclear 
%interactions that is absorbed locally
phi0=1000; %primary particle fluence
epsilon=0.2; %fraction of peak fluence in tail fluence

% testing functions
R=0:R0/100:R0*1.1; %cm
d=0:R0/100:R0*1.1; %cm
% R=0:0.1:R0;
% d=0:0.1:R0;

flu=fluence(phi0,beta,R0,d);
D_z=dose_C(phi0,beta,alpha,gamma,E0,p,d,rho,0,0,1);
Dhat_z=dosehat(phi0,beta,alpha,gamma,E0,p,d,rho,epsilon);

% % basic plots
% figure(3);
% plot(d,D_z./flu); hold on;
% ylabel('Dose per fluence');
% plot(d,Dhat_z./flu);
% figleg=legend('D(z)','$\hat{D}$(z)'); set(figleg,'Interpreter','Latex');
% figtitle=title(strcat('$\epsilon$ = ',num2str(epsilon))); set(figtitle,'Interpreter','Latex');


% % checking C library against Matlab function for consistency
% D=dose(phi0,beta,alpha,gamma,E0,p,d,rho,0,0,0); %non-C function
figure(4);
% plot(d,D./flu); hold on;
ax=plotyy(d,D_z./flu,d,D_z*rho./flu); hold on;
ylab=ylabel(ax(1),'Dose per fluence ($MeV g^{-1} cm^{2}$)'); set(ylab,'Interpreter','Latex');
ylab2=ylabel(ax(2),'LET ($MeV cm^{-1}$)');set(ylab2,'Interpreter','Latex');
xlabel('Depth (cm)');
figtitle=title(strcat('E = ',num2str(E0),'MeV')); set(figtitle,'Interpreter','Latex');
% figleg=legend('Matlab function','C Library'); set(figleg,'Interpreter','Latex');

% % % plotting effects of beam energy spread;
E_sigma=0.01*E0;
% D_realbeam=dose_C(phi0,beta,alpha,gamma,E0,p,d,rho,epsilon,E_sigma,1);
% figure(5);
% plot(d,D_realbeam./flu);hold on;
% plot(d,D_z./flu);
% figtitle=title(strcat('E = ',num2str(E0),'MeV')); set(figtitle,'Interpreter','Latex');
% figleg=legend(strcat('$\epsilon$ = ',num2str(epsilon),', $\sigma_{E,0}$ = ',num2str(E_sigma)),'Monochromatic beam'); set(figleg,'Interpreter','Latex');

% local heating
r1=3; %tube radius, in nm
r2=10;
D_heating=dose_C(phi0,beta,alpha,gamma,E0,p,d,rho,epsilon,E_sigma,1);
D_heating=D_heating./flu; %MeV cm^2 g^-1
LET=D_heating.*rho; %MeV cm^-1
dT1=(LET.*1e6*1.6e-19*1e-4)/(rho*pi*(r1*1e-9)^2*4.2); 
dT2=(LET.*1e6*1.6e-19*1e-4)/(rho*pi*(r2*1e-9)^2*4.2);
figure(6);
plot(d,dT1); hold on;
plot(d,dT2);
ylab=ylabel('$\Delta$ T/K'); set(ylab,'Interpreter','Latex');
xlabel('Depth/m'); 
legend(strcat('r = ',num2str(r1),'nm'), strcat('r =',num2str(r2),'nm'));
figtitle=title(strcat(num2str(E0),' MeV proton'));
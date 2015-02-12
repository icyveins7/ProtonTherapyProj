% Based on Bortfeld 1997: An analytical approximation of the Bragg curve
% for therapeutic proton beams
% Last edited: 22nd Jan
% Command line file to adjust parameters, plot graphs from functions

close all; clear all; clc; clear mem;

% Declaring constants/parameters
E0=100; %E0 in units of MeV
alpha=2.2e-3;
p=1.77;
rho=1; %mass density of medium, g/cm^3
R0=range(alpha,E0,p); %R in units of cm
beta=0.012;
gamma=0.6; %fraction of energy released in the nonelastic nuclear 
%interactions that is absorbed locally
phi0=1000; %primary particle fluence
epsilon=0.1; %fraction of peak fluence in tail fluence

% testing functions
R=0:R0/500:R0*1.1; %cm
d=0:R0/500:R0*1.1; %cm
% R=0:0.1:R0;
% d=0:0.1:R0;

flu=fluence(phi0,beta,R0,d);
D_z=dose_C(phi0,beta,alpha,gamma,E0,p,d,rho,0,0,1);
Dhat_z=dosehat(phi0,beta,alpha,gamma,E0,p,d,rho,epsilon);

% % basic plots
% figure(3);
% plot(d,D_z./flu); hold on;
% ylab=ylabel('Dose per fluence (MeV $cm^2$ $g^{-1}$)');set(ylab,'Interpreter','Latex');
% plot(d,Dhat_z./flu,'--');
% figleg=legend('D(z)','$\hat{D}$(z)'); set(figleg,'Interpreter','Latex');
% xlabel('Depth (cm)');
% % figtitle=title(strcat('$\epsilon$ = ',num2str(epsilon))); set(figtitle,'Interpreter','Latex');

% % energy comparisons
% figure(9);
% D_z2=dose_C(phi0,beta,alpha,gamma,0.75*E0,p,d,rho,0,0,1);
% D_z3=dose_C(phi0,beta,alpha,gamma,0.5*E0,p,d,rho,0,0,1);
% plot(d,D_z./flu); hold on;
% ylab=ylabel('Dose per fluence (MeV $cm^2$ $g^{-1}$)');set(ylab,'Interpreter','Latex');
% plot(d,D_z2./flu,'--'); plot(d,D_z3./flu,'k.-');
% figleg=legend(strcat(num2str(E0),' MeV'),strcat(num2str(0.75*E0),' MeV'),strcat(num2str(0.5*E0),' MeV')); set(figleg,'Interpreter','Latex');
% xlabel('Depth (cm)');
% % figtitle=title(strcat('$\epsilon$ = ',num2str(epsilon))); set(figtitle,'Interpreter','Latex');


% % % checking C library against Matlab function for consistency
% % D=dose(phi0,beta,alpha,gamma,E0,p,d,rho,0,0,0); %non-C function
% figure(4);
% % plot(d,D./flu); hold on;
% ax=plotyy(d,D_z./flu,d,D_z*rho./flu); hold on;
% ylab=ylabel(ax(1),'Dose per fluence ($MeV g^{-1} cm^{2}$)'); set(ylab,'Interpreter','Latex');
% ylab2=ylabel(ax(2),'LET ($MeV cm^{-1}$)');set(ylab2,'Interpreter','Latex');
% xlabel('Depth (cm)');
% figtitle=title(strcat('E = ',num2str(E0),'MeV')); set(figtitle,'Interpreter','Latex');
% % figleg=legend('Matlab function','C Library'); set(figleg,'Interpreter','Latex');

% % % plotting effects of beam energy spread;
E_sigma=0.01*E0;
% D_realbeam=dose_C(phi0,beta,alpha,gamma,E0,p,d,rho,epsilon,E_sigma,1);
% D_realbeam2=dose_C(phi0,beta,alpha,gamma,E0,p,d,rho,2*epsilon,2*E_sigma,1);
% figure(5);
% plot(d,D_realbeam./flu,'o','MarkerSize',1);hold on;
% plot(d,D_realbeam2./flu,'k-');
% plot(d,D_z./flu,'--');
% % figtitle=title(strcat('E = ',num2str(E0),'MeV')); set(figtitle,'Interpreter','Latex');
% figleg=legend(strcat('$\epsilon$ = ',num2str(epsilon),', $\sigma_{E,0}$ = ',num2str(E_sigma),' MeV'),strcat('$\epsilon$ = ',num2str(2*epsilon),', $\sigma_{E,0}$ = ',num2str(2*E_sigma),' MeV'),'Monochromatic beam','Location','NorthWest'); set(figleg,'Interpreter','Latex');
% ylab=ylabel('Dose per fluence (MeV $cm^2$ $g^{-1}$)');set(ylab,'Interpreter','Latex');
% xlabel('Depth (cm)');


% local heating
% r1=3; %tube radius, in nm
% r2=5;
% r3=7;
% r4=10;
% D_heating=dose_C(phi0,beta,alpha,gamma,E0,p,d,rho,epsilon,E_sigma,1);
% D_heating=D_heating./flu; %MeV cm^2 g^-1
% LET=D_heating.*rho; %MeV cm^-1
% dT1=(LET.*1e6*1.6e-19*1e-4)/(rho*pi*(r1*1e-9)^2*4.2); 
% dT2=(LET.*1e6*1.6e-19*1e-4)/(rho*pi*(r2*1e-9)^2*4.2);
% dT3=(LET.*1e6*1.6e-19*1e-4)/(rho*pi*(r3*1e-9)^2*4.2);
% dT4=(LET.*1e6*1.6e-19*1e-4)/(rho*pi*(r4*1e-9)^2*4.2);
% figure(6);
% plot(d,dT1); hold on;
% plot(d,dT2,'--'); plot(d,dT3,'o', 'MarkerSize',2); plot(d,dT4,'x', 'MarkerSize',2);
% ylab=ylabel('$\Delta$ T (K)'); set(ylab,'Interpreter','Latex');
% xlabel('Depth (cm)'); 
% legend(strcat('r = ',num2str(r1),'nm'), strcat('r =',num2str(r2),'nm'), strcat('r =',num2str(r3),'nm'), strcat('r =',num2str(r4),'nm'), 'Location','NorthWest');
% figtitle=title(strcat(num2str(E0),' MeV proton'));

% Gaussian noise 
num=1000;
E0values=normrnd(E0, 0.01*E0, [1,num]);
D_noise=zeros(num,length(d));
figure(7)
for i=1:num
    D_noise(i,1:end)=dose_C(phi0,beta,alpha,gamma,E0values(i),p,d,rho,0,0,1);
%     plot(d,D_noise(i,1:end)); hold on;
end
plot(d,sum(D_noise)/num,'k'); hold on;
dcompare=dose_C(phi0,beta,alpha,gamma,E0,p,d,rho,0,0.01*E0,1);
plot(d,dcompare);

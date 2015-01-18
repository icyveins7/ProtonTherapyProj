% Based on Bortfeld 1996: An analytical approximation of depth-dose distributions
% for therapeutic proton beams
% Last edited: 16th Jan
% Command line file to adjust parameters, plot graphs from functions

close all; clear all; clc; clear mem;

% declaring constants/parameters
E0=100; %E0 in units of MeV
alpha=2.2e-3;
p=1.77;
rho=1;
maxR=range(alpha,E0,p); %R in units of cm
d_b=maxR;
d_a=0.5*maxR;

% testing functions
R=0:0.01:maxR*1.1;
d=0:0.01:maxR*1.1;
D_BP_d=depthdose(alpha,E0,p,d,rho);
D0=max(D_BP_d);

figure(3);
plot(d,D_BP_d); hold on;
plot(d,abs(D_BP_d)); % probably not this one, according to paper
plot([maxR, maxR], get(gca, 'ylim'),'--');
legend('real','abs');

W_R=weight(alpha,p,R,D0,rho,d_a,d_b);
figure(4);
plot(R,W_R);hold on;
plot([maxR, maxR], get(gca, 'ylim'),'--');

W_disc_R=weight_disc(alpha,p,D0,rho,d_a,d_b,0.1);
figure(5);
plot(W_disc_R(2,1:end),W_disc_R(1,1:end),'.-'); hold on; %goes negative??
plot([maxR, maxR], get(gca, 'ylim'),'--');

% making discrete intervals of dose distribs
dim=size(W_disc_R);
dim=dim(2);
depthdose_mat=zeros(dim,length(d));
figure(6);
for i=1:dim
    E0value=(W_disc_R(2,i)/alpha)^(1/p);
    depthdose_mat(i,1:end)=depthdose(alpha,E0value,p,d,rho)*W_disc_R(1,i);
    plot(d,depthdose_mat(i,1:end)); hold on;
end
plot(d,sum(depthdose_mat,1));

% figure(7); % as a check for the raw dose distributions (unweighted)
% for i=1:dim
%     plot(


% % test convolution
% dsobp=conv(D_BP_d,W_R);
% figure(3);
% plot(d,dsobp);
% Based on Toulemonde, Surdutovich, Solov'yov (2009)
% Last edited by Ping Lin/Gabriel, 27th Jan
% Command line file to solve and plot temperature against time

close all; clear all; clc; clear mem;

currfolder=pwd;
currfolder=currfolder(1:end-15);
bortfolder=strcat(currfolder,'Bortfeld_1997');
addpath(bortfolder);

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

% Energy density tests
% r=logspace(-7,0,250); % mm
% [e1,e1_1]=energydensity_r(r,1);
% [e10,e10_1]=energydensity_r(r,10);
% [e20,e20_1]=energydensity_r(r,20);
% [e50,e50_1]=energydensity_r(r,50);
% [e100,e100_1]=energydensity_r(r,100);
% figure(3);
% loglog(r*1e6,e1,'b'); hold on; loglog(r*1e6,e1_1,'b--');
% loglog(r*1e6,e10,'r'); loglog(r*1e6,e10_1,'r--');
% loglog(r*1e6,e20,'Color',[0,0.5,1]); loglog(r*1e6,e20_1,'Color',[0,0.5,1],'LineStyle','--');
% loglog(r*1e6,e50,'k'); loglog(r*1e6,e50_1,'k--');
% loglog(r*1e6,e100,'m'); loglog(r*1e6,e100_1,'m--');
% legend('1 MeV','1 MeV(without correction)','10 MeV','10 MeV(without correction)','20 MeV','20 MeV(without correction)','50 MeV','50 MeV(without correction)','100 MeV','100MeV(without correction)');
% ylim([1e-8, 1e7]);

% // all code between // used in fcoefffunction, just typed here for
% // checking purposes
% solving for normalisation constant
% integrate exp(-(t-t0)^2/(2*s^2)) from 

% WANT TO WORK IN FOLLOWING UNITS, convert anything to these
% J g K nm s

E=10; %MeV
fun=@(r)energydensity_r(r,E).*r; % r in mm
r_integral=integral(fun,0,1); % J kg^-1 mm^2
r_integral=r_integral*1e-3*1e12; % J g^-1 nm^2
t_integral=integral(@energydensity_t,0,1e-9); % integrate from - inf instead??
bconst=1/(t_integral*r_integral*2*pi);
% need to solve for value of E of proton at the bragg peak->solve for v ->
% insert into energydensity_r -> use LET derived from Bortfeld funcs
% copied
alpha=2.2e-3;
p=1.77;
rho=1; %mass density of medium, g/cm^3
R0=range(alpha,E,p); %R in units of cm
beta=0.012;
gamma=0.6; %fraction of energy released in the nonelastic nuclear 
%interactions that is absorbed locally
phi0=1000; %primary particle fluence

depth=0:R0/100:R0*1.1; %cm

% code to solve remaining proton energy at bragg peak, pasted into
% fcoefffunction
D=dose_C(phi0,beta,alpha,gamma,E,p,depth,rho,0,0,1); % MeV/g
% find LET of bragg peak
flu=fluence(phi0,beta,R0,depth);
LET=D*rho./flu; % MeV/cm, for water it's the same as dose per fluence since rho=1g/cm^3
LET = LET * 1.602e-19 * 1e6/1e7; % J/nm
[S_e,ind]=max(LET); % checked against SRIM, correct
figure;
ax=plotyy(depth,D,depth,LET);
ylab=ylabel(ax(1),'Dose per fluence ($MeV g^{-1} cm^{2}$)'); set(ylab,'Interpreter','Latex');
ylab2=ylabel(ax(2),'LET ($J nm^{-1}$)');set(ylab2,'Interpreter','Latex');

% find remaining energy of proton at bragg peak, this is used as input to
% energydensity_r later, to solve for velocity or beta
rem_E=E*1e6*1.602e-19 -trapz(depth(1:ind).*1e7,LET(1:ind)); % J
rem_E_MeV=rem_E/(1.602e-19)/1e6;
% //

% override values
S_e=0.08*1e3*1.602e-19;
rem_E_MeV=0.04;

% % calculate temperature spike at bragg peak
% radius=25;
% % make geometry description matrix
% gd=[1,0,0,radius]'; % circle with specified radius (in nm) centred at 0,0
% g=decsg(gd,'C1',('C1')');
% % make mesh
% [p,e,t] = initmesh(g,'hmax', radius,'MesherVersion','R2013a');
% for i=1:6
%     [p,e,t] = refinemesh(g, p, e, t);
% end
% % check plots
% figure;
% pdegplot(g,'edgeLabels','on');
% figure;
% pdeplot(p,e,t);
% % Create a pde entity for a PDE with 2 dependent variables
% numberOfPDE = 2;
% pb = pde(numberOfPDE);
% % Create a geometry entity
% pg = pdeGeometryFromEdges(g);
% bOuter = pdeBoundaryConditions(pg.Edges(1:4), 'u', [310,310]'); % 310K boundary conditions, both systems
% pb.BoundaryConditions = bOuter;
% 
% tlist=logspace(-16,-9);
% c=[2e-7;0;2e-7;6e-10;0;6e-10]; % [Ke,K] in units of W K^-1 nm^-1 
% lmbda=2; % nm
% g_const=2e-7/(lmbda^2); % e-phonon coupling
% a= [g_const, -g_const, -g_const, g_const]';
% % Ce, in J K^-1 g^-1, rho*C, in g nm^-3 J g^-1 K^-1
% d=[3.0468e-25 0 0 4.18*rho*1e-21]'; 
% % init conditions
% % Assume N and p exist; N = 1 for a scalar problem
% N=2;
% np = size(p,2); % Number of mesh points
% u0 = zeros(np,N); % Allocate initial matrix
% for k = 1:np
%     x = p(1,k);
%     y = p(2,k);
%     u0(k,:) = 310; % Fill in row k, constant temperature throughout, 310K
% end
% u0 = u0(:); % Convert to column form
% f_placeholder=@(p,t,u,time)fcoeffunction(p,t,u,time,rem_E_MeV,bconst,S_e);
% u = parabolic(u0,tlist,pb,p,e,t,c,a,f_placeholder,d);
% % plot the molecular system 
% figure;
% pdeplot(p,e,t,'xydata',u(np*N/2+1:end,end),'contour','on','colormap','hot');
% title('Molecular system, t=end');
% figure;
% pdeplot(p,e,t,'xydata',u(np*N/2+1:end,26),'contour','on','colormap','hot');
% title('Molecular system, t~1e-14s');
% % % plot the electronic system
% % figure;
% % pdeplot(p,e,t,'xydata',u(1:np*N/2,end),'contour','on','colormap','hot');
% % title('Electronic system');

% % molecular system temp graphs
% % find middle point, 1 nm point etc.
% rlist=(p(1,:).^2+p(2,:).^2).^0.5;
% figure;
% legentries={};
% for i=0:0.5:5
%     [nil,ind]=min(abs(rlist-i));
%     semilogx(tlist,u(np+ind,:)); hold on;
%     legentries{end+1}=strcat(num2str(rlist(ind)),' nm');
% end
% xlabel('Time (s)'); ylabel('Temp. (K)'); 
% legend(legentries);

lmbda=2; % nm
g_const=2e-7/(lmbda^2); % e-phonon coupling
% g_const=g_const/25;
pdefun=@(x,t,u,dudx)pdefun_plchold(x,t,u,dudx,rho,g_const,rem_E_MeV,bconst,S_e);
icfun=@(x) [310;310]; % initial conditions 310K both systems
space_itv=0.1;
xmesh=0:space_itv:25;
tspan=logspace(-16,-8,250);
sol=pdepe(1,pdefun,icfun,@bcfun,xmesh,tspan);
figure;
legentries={};
for i=0:0.5:5
    [nil,ind]=min(abs(xmesh-i));
    semilogx(tspan,sol(:,ind,2)); hold on;
    legentries{end+1}=strcat(num2str(xmesh(ind)),' nm');
end
xlabel('Time (s)'); ylabel('Temp. (K)'); 
legend(legentries); grid on;

% % electronic system
% figure;
% legentries_e={};
% for i=0:0.5:5
%     [nil,ind]=min(abs(xmesh-i));
%     semilogx(tspan,sol(:,ind,1)); hold on;
%     legentries_e{end+1}=strcat(num2str(xmesh(ind)),' nm');
% end
% xlabel('Time (s)'); ylabel('Temp. (K)'); 
% legend(legentries_e); grid on;


% % make contour from radial pdepe solution
% % start by making radial mesh, here we use m*n points per circle for the
% % nth circle from centre
% m=8; t_ind=108;
% % find max temperature of that time to define colours
% maxtemp=max(sol(t_ind,:,2));
% figure(5);
% for i=length(xmesh):-1:1
%     XCON=[]; YCON=[]; TCON=[];
%     for j=1:m*(i-1)
%         angle=(j/(m*(i-1)))*2*pi;
%         XCON(end+1)=xmesh(i)*cos(angle);
%         YCON(end+1)=xmesh(i)*sin(angle);
%         TCON(end+1)=sol(t_ind,i,2);
%     end
%     if i~=1
%         %define contour colours, all in this contour should be same temperature
%         if TCON(1)<310+(maxtemp-310)/3
%             colour=[(TCON(1)-310)/((maxtemp-310)/3), 0, 0];
%         elseif TCON(1)<310+2*(maxtemp-310)/3
%             colour=[ 1, (TCON(1)-(310+(maxtemp-310)/3))/((maxtemp-310)/3), 0];
%         else
%             colour=[ 1, 1, (TCON(1)-(310+2*(maxtemp-310)/3))/((maxtemp-310)/3)];
%         end
%         fill(XCON,YCON,TCON,'FaceColor',colour,'EdgeColor',colour); hold on;
%     end
% end
% title(strcat('Time: ',num2str(tspan(t_ind)),' s'));

% % make movie from frames
% m=8; frames_array(length(tspan))=struct('cdata',[],'colormap',[]);
% fig=figure('Position',[100,100,1000,900]); %setting resolution
% for t_ind=1:length(tspan)
%     % find max temperature of that time to define colours
%     % maxtemp=max(sol(t_ind,:,2));
%     maxtemp=550;
%     
%     for i=length(xmesh):-1:1
%         XCON=[]; YCON=[]; TCON=[];
%         for j=1:m*(i-1)
%             angle=(j/(m*(i-1)))*2*pi;
%             XCON(end+1)=xmesh(i)*cos(angle);
%             YCON(end+1)=xmesh(i)*sin(angle);
%             TCON(end+1)=sol(t_ind,i,2);
%         end
%         if i~=1
%             %define contour colours, all in this contour should be same temperature
%             if TCON(1)<=310
%                 colour=[0,0,0];
%             elseif TCON(1)<310+(maxtemp-310)/3
%                 colour=[(TCON(1)-310)/((maxtemp-310)/3), 0, 0];
%             elseif TCON(1)<310+2*(maxtemp-310)/3
%                 colour=[ 1, (TCON(1)-(310+(maxtemp-310)/3))/((maxtemp-310)/3), 0];
%             else
%                 colour=[ 1, 1, (TCON(1)-(310+2*(maxtemp-310)/3))/((maxtemp-310)/3)];
%             end
%             fill(XCON,YCON,TCON,'FaceColor',colour,'EdgeColor',colour); hold on;
%         end
%     end
%     title(strcat('Time: ',num2str(tspan(t_ind)),' s'));
%     frames_array(t_ind)=getframe(fig);
%     clf(fig);
% end
% % save movie
% save('movie.mat','frames_array');

% load movie
load('movie.mat');
% play the movie
mov_fig=figure('Position',[100,100,1000,900]); %setting resolution
movie(mov_fig,frames_array,5,20)

rmpath(bortfolder);

% testing
E_ion=1; %MeV
lorentzfactor=(E_ion/938)+1;
b=(1-(1/lorentzfactor)^2)^0.5; %beta
capW=(2*9.10938291e-31*(299792458)^2*b^2*(1-b^2)^-1)/(1.60217657e-16); % kinematically limited max energy, in keV
k=6e-11;
density=1e-6;
r=1e-6;
lowerlimit=(r*density/k)^(1/1.079);
% lowerlimit=(r*density/k)^(1/1.667);
answer=integral(@testintegral,lowerlimit,capW,'RelTol',1e-15,'AbsTol',1e-15);
wlist=lowerlimit:(capW-lowerlimit)/100000:capW;
plot(wlist,testintegral(wlist)); hold on;
% plot(wlist,imag(testintegral(wlist)));
answer_trapz=trapz(wlist,testintegral(wlist));

% for alpha=1.079:0.001:1.667
%     lowerlimit=(1e-6*density/k)^(1/alpha);
%     answer=integral(@testintegral,lowerlimit,capW,'RelTol',1e-9,'AbsTol',1e-12);
%     plot(alpha,answer,'x'); hold on;
%     plot(alpha,[9.9727e11],'s');
% end

Rmax= k*capW^1.667;
answer_real=((1-(r*density/Rmax))^(1/1.079)) / (r*density);

r0=k*0.078^1.079;
answer_realwithI= (1-(r*density + r0)/(Rmax + r0))^(1/1.079) / (r*density + r0);
% answer_quad=quad(@testintegral,lowerlimit,capW);

c1=1.352817016; % Ne^4/mc^2 in units of keV mm^-1
c1=c1*(1.60217657e-16); % J mm^-1

% constants for Rudd, arranged in ascending ionisation potential order
A1=[1.02 1.02 1.02 1.02 1.25]; A2=[1.07 1.07 1.07 1.07 1.1];
B1=[82 82 82 82 0.5]; B2=[14.6 14.6 14.6 14.6 1.3];
C1=[0.45 0.45 0.45 0.45 1]; C2=[0.6 0.6 0.6 0.6 1];
D1=[-0.8 -0.8 -0.8 -0.8 1]; D2=[0.04 0.04 0.04 0.04 0];
E1=[0.38 0.38 0.38 0.38 3]; alpha=[0.64 0.64 0.64 0.64 0.66];
G=[0.99 1.11 1.11 0.52 1];

I=[12.61 14.73 18.55 32.2 539.7]; %Binding energy of the shell (water vapour) (eV)
Il=[10.79 13.39 16.05 32.3 539]; % Binding energy for liquid water (eV)

a0=0.0529*10^-9;  %Bohr radius (m)
R=13.606;     %Rydberg constant (eV)
N=2;        %Shell occupancy
m=9.109*(10^-31); %Electron's mass (kg)
M=(1.673*10^-27); %Proton's mass (kg)
w=[W/I(1); W/I(2); W/I(3); W/I(4); W/I(5)]; %Dimensionless normalized kinetic energy of the ejected electron
beta=sqrt(1-(M*(2.998*10^8)^2/(M*(2.998*10^8)^2+E_ion*1.602*10^-19))^2); % Relativistic effects
c=299792458;
z=1;

for i=1:5
    v=sqrt(m*(beta*c)^2/(2*I(i)*1.602e-19)); %Dimensionless normalized velocity
    F1=A1(i)*(log((1+v(i)^2)/(1-beta^2))-beta^2)/(B1(i)/v(i)^2+v(i)^2)+(C1(i)*v(i)^D1(i))/(1+E1(i)*v(i)^(D1(i)+4));        
    F2=C2(i)*(v(i)^D2(i))*(A2(i)*v(i)^2+B2(i))/(C2(i)*v(i)^(D2(i)+4)+A2(i)*v(i)^2+B2(i));
    wmax(i)=4*(v(i)^2)-2*v(i)-R/(4*I(i));
    tcs =G(i).*((z^2)*4*pi*(a0^2)*N*(R^2)/I(i)^3).*((F1+w(i)*F2)./(((1+w(i)).^3).*(1+exp(alpha(i).*(w(i)-wmax(i))./v(i)))))'; %sdcs
end
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

% % Energy density tests
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
function [x,L_E_T]=BPprotons(Tn, namefile, inc_analytic)
% Adding default arguments to old code; code should run as before with 2
% arguments, but inc_analytic can be included to plot the Bortfeld_1997
% approximations if set to 1.
% Edited 25 Jan 2015 by Gabriel/Pinglin
switch nargin
    case 2
        inc_analytic=0;
end
%we we want to run the program we write for example that Bpprotons(100,'aslanin.txt') 
% and 100 means 100 MeV and aslanin.txt is the name of the file which
% will be made and its output will be read by SOBP code later.
% One feaure of our current code is that it also uses the matlab figure which
% has experimental data and then it is plotted here for comparison
% purposes.

% Tn= Initial kinetic energy of the proton, MeV
% L_E_T= Vector with the Linear Energy Transfer after considering straggling and fragmentation
% x = Vector with depths

% Range of energies for which x and L_E_T is calculated
T=[Tn:-0.5:50 49:-0.1:5 4.8:-0.2:3.2 3.1:-0.05:1.1 1:-0.01:0.01 0.009:-0.001:0.001 0.0009:-0.0001:0.0001]*10^6; dimT=length(T)

%According to dingfelder, Eq.57, for Stopping cross section we should consider 4 facotrs:
%1-ionization by proton,2-excitation by proton,3-ionizatio by hydrogen
%4-charge changing( i.e. charge transfer(electron capture S10)  and   stripping(electron loss S01)  altogether ).
LET=zeros(1,dimT); x=zeros(1,dimT); Stopcs=zeros(1,dimT); Stcshyd=zeros(1,dimT); StopCC=zeros(1,dimT); Stexc=zeros(1,dimT); % Wave=zeros(1,dimT);
Tcs=zeros(1,dimT); Tcshyd=zeros(1,dimT); S01=zeros(1,dimT); S10=zeros(1,dimT); shelltot=zeros(1,dimT);

% Here we call function let which calculates linear energy transfer but, 
% without straggling and fragmentation. Of course the LET is calcualted
% for the T(1) which is in fact initial energy of the projectile.
% The LET calculated here later will be put into another function L_E_T
% which adds straggling and fragmentation to the LET, calculated here and
% then this final LET as the Bragg peak is plotted.


% we call in the following the function "let" for T(1) and later we call "let" for T(k)
% for k from 2 until the end.

[LET(1),Stopcs(1),Stcshyd(1),Stexc(1),StopCC(1),Tcs(1),Tcshyd(1),S01(1),S10(1),shelltot(1)]=let(T(1));
for k=2:dimT
    % just to see if our programm is running well:
    if k==100 || k==200 || k==300 || k==400 || k==500 || k==600
        k
    end
    
     %we call "let" for T(k) for k from 2 until the end.
   [LET(k),Stopcs(k),Stcshyd(k),Stexc(k),StopCC(k),Tcs(k),Tcshyd(k),S01(k),S10(k),shelltot(k)]=let(T(k)); %Linear energy transfer

    
    % LET is calculated above(of course without straggling and
    % fragmentation)and now the positions should be calculated according to equation 10 of Obolensky paper.
    
    % Here in order to speed up the calculation of the integral, the total intergral
    % is calcualted between two successive steps and then all results are
    % added to each other. x(T)=int^{T0}_{T}(dT'/dT'/dx)

    Tp=[T(k) T(k-1)]; LETp=[LET(k) LET(k-1)]; y=1./LETp;
    x(k)=trapz(Tp,y)+x(k-1); %Depth m
end

figure;
loglog(T,Tcs, T,Tcshyd, T, S10, T, S01, T, shelltot);
xlabel('incidental particle energy')
ylabel('total cross section')

x=x*10^3; %mm   x in unit of meter is converted to x in units of cm.
LET=LET*10^-9; % MeV/mm

% % Plot Stopping cross section
% 
% % a matlab figure containing experimental result for total Stopping cross
% % section now opens
% % to be comapred with our own results.
% open ('Exper_stop_cross_section_prot_liq_water.fig');
% hold on; semilogx(T,Stopcs+Stexc+Stcshyd+StopCC,T,Stopcs,T,StopCC,T,Stcshyd,T,Stexc)
% 
% xlabel('Kinetic energy of the projectile, T (eV)'); ylabel('Stopping cross sections (eV m^2)');
% title('\it{Stopping cross sections for liquid water}','FontSize',16); 
% legend('Experimental','Total','Ionisation (p)','Charge Changing (p)','Ionisation (H)','Excitation (p)',1);


% Some extra points are added at the end of the peak with LET=0, with
% physical sense as we consider up to here that no particles will reach a
% depth larger than the Range, where the peak is located.
x=[x x(dimT)+0.001 x(dimT)+0.01 x(dimT)+0.05 x(dimT)+0.1 x(dimT)+0.2 x(dimT)+0.5 x(dimT)+1 x(dimT)+1.3 x(dimT)+1.6 x(dimT)+2 x(dimT)+2.3 x(dimT)+2.6 x(dimT)+3 x(dimT)+3.5 x(dimT)+4 x(dimT)+4.5 x(dimT)+5 x(dimT)+5.5 x(dimT)+6 x(dimT)+6.5 x(dimT)+7 x(dimT)+7.5 x(dimT)+8];
%the number of points we have added our length(x)-dim(T)=(length of the current x minus length of the former x)
% and we should consdier them to be zero.
LET=[LET zeros(1,length(x)-dimT)];

% Now the LET which was calculated prevoiusly without straggling and
% fragmentation is put into the following function L_E_T  and using Kundrat paper 
%straggling and fragmentatin are added and then the final
% LET which is what we want can be plotted for the Bragg peak. Note that we
% use particle numbers and not the fluence as Kundrat suses, therfore there
% is no need to be worried abut the wrong dimensions!
L_E_T=straggling(x,LET); % MeV/mm

% figure(1);
% plot(x,LET); xlabel('Depth, x (mm)'); ylabel('LET (MeV/mm)');

% Addition here by Gabriel/Pinglin
% Argument inc_analytic determines whether to call Bortfeld functions
% Change addpath and rmpath to folder location as necessary
figure(2); plot(x,L_E_T); xlabel('Depth, x (mm)'); ylabel('LET (MeV/mm)'); hold on;
if (inc_analytic~=0)
    disp('Calling Bortfeld_1997 analytic function...');
    addpath I:/Work/ProtonTherapyProj/Bortfeld_1997
    % we use default values for water, imported from Bortfeld commandline code
    % here x is in mm but Bortfeld code is in cm
    alpha=2.2e-3;
    p=1.77;
    rho=1; %mass density of medium, g/cm^3
    R0=range(alpha,Tn,p); %R in units of cm
    beta=0.012;
    sigma=0.012*(R0^0.951); % changed power to 0.951, same as this code
    gamma=0.6; %fraction of energy released in the nonelastic nuclear 
    %interactions that is absorbed locally
    phi0=1; %primary particle fluence
    epsilon=0; %fraction of peak fluence in tail fluence
    % calculate
    L_E_Tbf=dose_C(phi0,sigma,beta,alpha,gamma,Tn,p,x/10,rho,epsilon,1);
    L_E_Tbf=L_E_Tbf./(fluence(phi0,beta,R0,x/10)/100); %MeV mm^2 g^-1
    L_E_Tbf=L_E_Tbf.*rho*1e-3; %MeV mm^-1
    % plot
    plot(x,L_E_Tbf); legend('Candela code','Bortfeld code');
    rmpath I:/Work/ProtonTherapyProj/Bortfeld_1997
end
% end of additions

%SAVE x and L_E_T in a txt-file.
myfile = fopen(namefile,'w');
y=[x; L_E_T];
fprintf(myfile,'x (mm) \t LET (MeV/mm) \n');
fprintf(myfile,'%4.4E %4.4E\n', y);
fclose(myfile);

%------------------------------------------------------------------------
%_____In fact our  program is now finished_______
% and in the rest we just write the functions
%________________________________________________
%__________________________________________________________________________
%__________________________________________________________________________

%-----------------------------------------------------------------------
function L_E_T=straggling(x,LET)

% Function which adds STRAGGLING and Nuclear Fragmentation given a vector
% with the depths and their corresponding LET using Continuous Slowing Down
% Approximation (CSDA)
A=1;
N0=1; % Number of Incident particles
lambda=435; %mm Mean Free Path Length for nuclear reactions
      
[LETmax,index]=max(LET); %Determine the position and the value of the maximum of the LET.
R0=x(index); %mm  maximum penetration depth of the projectile
Sigma_Str=0.012*((R0/10)^0.951)*A^(-0.5)*10; %mm
%Sigma_Str=0.012*((x/10).^0.951)*A^(-0.5)*10; %mm

% x from the start of the program is in unit of mm, so is R0. 
% R0 divided by 10 cahnges the units of R0 from mm to cm as R0 should be in cm for
% calculating the
% the Sigma according to Kundrat.
% We later multiply the whole Sigma by 10 at the end to change Sigma from cm to mm.
%Sigma_Str=0.012*((x/10).^0.951)*A^(-0.5)*10; %mm

L_E_T=zeros(1,length(x));
for k=1:length(x)
    if x(k)<=5 %Threshold considered below which there's no straggling
        L_E_T(k)=LET(k);
    else
        
        % LET without fragmentation and straggling is put into the
        % following formula(Equatios:1,5 and 6 of Kundrat)so that straggling and
        % fragmentation is added
        integral=trapz(x, LET.*exp(-(x-x(k)).^2./(2*Sigma_Str^2)));
        L_E_T(k)=N0*exp(-x(k)/lambda)*integral/(sqrt(2*pi)*Sigma_Str);
    end
end
%


%------------------------------------------------------------------------

function [LET,Stopcs,Stcshyd,Stexc,StopCC,Tcs,Tcshyd,S01,S10,shelltot]=let(T)
% This function calculates the stopping cross sections for the different processes considered and the LET.

% Constants
z=1;        %Projectile's charge
Z=10; %Number of electrons in the target material
a0=0.0529*10^-9;  %Bohr radius (m)
R=13.606;     %Rydberg constant (eV)
N=2;        %Shell occupancy
m=9.109*(10^-31); %Electron's mass (kg)
M=(1.673*10^-27); %Proton's mass (kg)
n=3.343*10^28; % Number density of molecules for liquid water (molecules/m^3)
I=[12.61 14.73 18.55 32.2 539.7]; %Binding energy of the shell (water vapour) (eV)
Il=[10.79 13.39 16.05 32.3 539]; % Binding energy for liquid water (eV)

% Experimental parameters
A1=[1.02 1.02 1.02 1.02 1.25]; A2=[1.07 1.07 1.07 1.07 1.1];
B1=[82 82 82 82 0.5]; B2=[14.6 14.6 14.6 14.6 1.3];
C1=[0.45 0.45 0.45 0.45 1]; C2=[0.6 0.6 0.6 0.6 1];
D1=[-0.8 -0.8 -0.8 -0.8 1]; D2=[0.04 0.04 0.04 0.04 0];
E1=[0.38 0.38 0.38 0.38 3]; alpha=[0.64 0.64 0.64 0.64 0.66];
G=[0.99 1.11 1.11 0.52 1]; %Partitioning factor to adjust the contributions of the different subshells to the results obtained from Born approx.



 
%beta=sqrt(1-(M*(2.998*10^8)^2/(M*(2.998*10^8)^2+T*1.602*10^-19))^2); % Relativistic effects
%gamma=1/sqrt(1-beta^2);
%zeff=z*(1-exp(-125*beta*(z^-(2/3))));




%IONISATION cross-section by protons and neutral hydrogen.
v=zeros(1,5); wmax=zeros(1,5); Stopcs=0; Stcshyd=0; %Tcs=0; Tcshyd=0;
g=0.8/((1+exp((log10(T)-4.2)/0.5)))+0.9; %For the ionisation by neutral hydrogen

NE=10000; % Number of W
E=zeros(1,NE);

Tcs=0; Integral_Wave=0; Tcshyd=0;
  for i=1:5  % Sum for the different sub-shells.
    if T>=Il(i)  %If the projectile has energy enough to produce ionisation...
        
        E(NE)=T; E(1)=Il(i);
        for cont=1:NE-2
            E(cont+1)=E(1)*exp(cont*log(E(NE)/E(1))/(NE-1)); %Kinetic energy of the ejected electrons
        end
        W(1,:)=(E(:)-Il(i)); 
        % Although we are working with liquid water parameters, but in the following 
       % we use "I" as defined for water vapor(and not Il), because in Rudd's
       % semi-imperical method ,data have been fitted based on water vapor and
       % nobody can do anything now. Of course the other parameters are for
       % liquid water. Refere to Dingfelder paper, page 267, second column
       % for more explanations.
        w=[W/I(1); W/I(2); W/I(3); W/I(4); W/I(5)]; %Dimensionless normalized kinetic energy of the ejected electron
        v(i)=sqrt(m*(T)/(M*I(i))); %Dimensionless normalized velocity
        wmax(i)=4*(v(i)^2)-2*v(i)-R/(4*I(i));
        F1=A1(i)*log(1+v(i)^2)/(B1(i)/v(i)^2+v(i)^2)+(C1(i)*v(i)^D1(i))/(1+E1(i)*v(i)^(D1(i)+4));
        
       % if we consider the F1 given by the Surdutovich paper, which has
       % some corrections compared to the Obolensky paper:
        %F1=A1(i)*(log((1+v(i)^2)/(1-beta^2))-beta^2)/(B1(i)/v(i)^2+v(i)^2)+(C1(i)*v(i)^D1(i))/(1+E1(i)*v(i)^(D1(i)+4));

        
        F2=C2(i)*(v(i)^D2(i))*(A2(i)*v(i)^2+B2(i))/(C2(i)*v(i)^(D2(i)+4)+A2(i)*v(i)^2+B2(i));
        tcs =G(i).*((z^2)*4*pi*(a0^2)*N*(R^2)/I(i)^3).*((F1+w(i,:)*F2)./(((1+w(i,:)).^3).*(1+exp(alpha(i).*(w(i,:)-wmax(i))./v(i)))))'; %sdcs
        tcshyd=g.*tcs; %sdcs for hydrogen
        
        %Sdcs=Sdcs+tcs'; Sdcshyd=Sdcshyd+tcshyd'; %Accumulative single differentiation cross section for protons and hydrogen
        
        Tcs=Tcs+trapz(E,tcs); %Total ionisation cross section (proton)
        Integral_Wave=Integral_Wave+trapz(E,E.*tcs'); %Integral to compute the W average
        Tcshyd=Tcshyd+trapz(W,tcshyd); %Total ionisation cross section (hydrogen)
        
        Stopcs=Stopcs+trapz(E,E.*tcs'); %Ionisation by protons
        Stcshyd=Stcshyd+trapz(E,E.*tcshyd'); %Ionisation by hydrogen molecule
    end
  end
 
Wave=Integral_Wave/Tcs; % Average energy of the secondary electrons


% EXCITATION cross section (by protons), Sexc
  % - Parameters (Dingfelder, Inokuti, Herwig, Paretzke. Inelastic-collision cross sections of liquid water for protons (p.266)
  
  %Here, we first calcualte the cross section for the  excitation  based on equation 34 of Dingfelder and call it
   % "shell" and then we calculate the stopping cross section for
   % excitation based on equation  15 of the second paper by Obolenski and
   % call it Stexc.
  Ek=[8.17 10.13 11.31 12.91 14.5]; sigma0=10^-20; %Parameters of lenght in m and parameters of energy in eV
  a=[876 2084 1373 692 900]; J=[19820 23490 27770 30830 33080];
  omega=[0.85 0.88 0.88 0.78 0.78]; nu=[1 1 1 1 1];
  %C1=[3.132 3.468 3.324 3.028 3.028]; C2=[0.25 0.25 0.25 0.25 0.25];
  %Tmax=[47.15 54.78 65.75 77.22 82.86]; sigmax=[8 1.3 8.11 5.19 6.03]*10^-22;
shelltot=0; Stexc=0;% Sexc=0;
for i=1:5
    shell=(sigma0*(Z*a(i))^omega(i)*(T-Ek(i))^nu(i))/(J(i)^(omega(i)+nu(i))+T^(omega(i)+nu(i)));
    %Sexc=Sexc+shell;
    %Tpri=0:0.01:Ek(i); %Sexcpri=zeros(1,length(Tpri));
    %for j=1:length(Tpri)
    %    Sexcpri(j)=(sigma0*(Z*a(i))^omega(i)*(Tpri(j)-Ek(i))^nu(i))/(J(i)^(omega(i)+nu(i))+Tpri(j)^(omega(i)+nu(i)));
    %end
     shelltot=shelltot+shell;
    Stexc=Stexc+shell*Ek(i);%-trapz(Tpri,Sexcpri);
end
% CHARGE TRANSFER cross section, S10
a00=-0.18; b0=-18.22; c0=0.215; d0=3.55; a1=-3.6; b1=-1.997; x0=3.45; x1=5.251; % Constant parameters to fit experimental values
X=log10(T);
Y=(a00*X+b0-c0*(X-x0)^d0*heaviside(X-x0))*heaviside(x1-X)+(a1*X+b1)*heaviside(X-x1);
S10=10^Y;

% STRIPPING cross section (of neutral hydrogen), S01
A=2.835; B=0.31; C=2.1; D=0.76; % Constant parameters to fit experimental values
tau=m*T/M; % Kinetic energy of an electron travelling with the same speed as the proton
slow=4*pi*(a0^2)*C*(tau/R)^D;
shigh=4*pi*(a0^2)*(R/tau)*(A*log(1+tau/R)+B);
S01=(slow*shigh)/(shigh+slow);

% PROBABILITY to find the ion in the charge state: neutral hydrogen probH or proton probp
probH=S10/(S01+S10);
probp=S01/(S01+S10);

%Total stopping cross section TCS
StopCC=((S10*S01)/(S01+S10))*(Il(1)+tau); %Charge-changing stopping cross section
Stopcs=probp*Stopcs; % Ionisation stopping cross section by protons      
Stcshyd=probH*Stcshyd; % Ionisation stopping cross section by Hydrogen   
Stexc=probp*Stexc;% Excitation cross section should also be given a % probability.

TCS=StopCC+Stopcs+Stcshyd+Stexc;
%Linear energy transfer (eV/m)
LET=n*TCS;

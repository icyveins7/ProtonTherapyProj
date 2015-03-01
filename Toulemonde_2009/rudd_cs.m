function value = rudd_cs(W,i,E_ion)

% units of differential cross section from old code appear to be SI, except
% for energy which is in eV i.e. final units in m^2/eV? we want it in
% mm^2/J to convert to Gy. older rutherford energydensity_r had it in units
% of mm and keV. so we want W in eV here, scale all conditions accordingly.
% E_ion accepted in MeV.

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
beta=sqrt(1-(M*(2.998*10^8)^2/(M*(2.998*10^8)^2+E_ion*1e6*1.602*10^-19))^2); % Relativistic effects
c=299792458;
z=1;


v=sqrt(m*(beta*c)^2/(2*I(i)*1.602e-19)); %Dimensionless normalized velocity
F1=A1(i)*(log((1+v^2)/(1-beta^2))-beta^2)/(B1(i)/v^2+v^2)+(C1(i)*v^D1(i))/(1+E1(i)*v^(D1(i)+4));        
F2=C2(i)*(v^D2(i))*(A2(i)*v^2+B2(i))/(C2(i)*v^(D2(i)+4)+A2(i)*v^2+B2(i));
wmax=4*(v^2)-2*v-R/(4*I(i));
value =G(i).*((z^2)*4*pi*(a0^2)*N*(R^2)/I(i)^3).*((F1+w(i,:).*F2)./(((1+w(i,:)).^3).*(1+exp(alpha(i).*(w(i,:)-wmax)./v)))); %sdcs


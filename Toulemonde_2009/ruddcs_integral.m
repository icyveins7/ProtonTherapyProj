function value=ruddcs_integral(W,r,i,E_ion)
% W in eV
% tinkering with rudd
% use r in m, NOT mm
n=3.343*10^28; % Number density of molecules for liquid water (molecules/m^3)
Ne=n*10; % number density of electrons (10 per molecule) for liquid water (electrons/m^3)

density=1e3; % kg m^-3
r=r*density; % m -> kg/m^2

alpha_w=zeros(1,length(W));

for j=1:length(W)
    if (W(j)<1000)
        alpha_w(j)=1.079;
    else 
        alpha_w(j)=1.667;
    end
end
k=6e-11*1e6; % g cm^-2 keV^-alpha_w -> kg mm^-2 keV^-alpha_w -> kg m^-2 keV^-alpha_w

% want to convert all to kg, m, eV, convert to J at the end, after
% integration
% rudd_cs accepts in eV, outputs cross section in m^2 eV^-1

value= (Ne/5).*rudd_cs(W,i,E_ion).* (W./(alpha_w.*k.*(W./1000).^(alpha_w))) .* (1 - r./(k.*(W./1000).^alpha_w)).^(1./alpha_w - 1);
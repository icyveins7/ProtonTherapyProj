function value=ruddcs_integral(W,r,i,E_ion)
% W in keV
% tinkering with rudd
n=3.343*10^28; % Number density of molecules for liquid water (molecules/m^3)
Ne=n*10/(1e9); % number density of electrons (10 per molecule) for liquid water (electrons/mm^3)

density=1e-6; % kg mm^-3
r=r*density; % mm -> kg/mm^2

if (W<1)
    alpha_w=1.079;
else 
    alpha_w=1.667;
end
k=6e-11; % g cm^-2 keV^-alpha_w -> kg mm^-2 keV^-alpha_w

% want to convert all to kg, mm, keV, convert to J at the end, after
% integration
% rudd_cs accepts in eV, outputs cross section in m^2 eV^-1, so convert to
% mm^2 keV^-1 -> *1e-6*1e3=1e-3
% we add the 1/(2pi*r) factor here which was included outside previously
% r are in units of mass/length^2, so total units here should be
% correct now i.e. energy/mass when integrated with dw, keV/kg
value= Ne.*rudd_cs(W*1000,i,E_ion).*1e-3.* (1./(alpha_w.*k.*W.^(alpha_w-1))) .* (1 - r./(k.*W.^alpha_w)).^(1./alpha_w - 1);
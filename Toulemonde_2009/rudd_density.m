function energydensity=rudd_density(rlist,E_ion)
% r in m, E_ion in MeV

lorentzfactor=(E_ion/938)+1;
b=(1-(1/lorentzfactor)^2)^0.5; %beta
capW=1e3*(2*9.10938291e-31*(299792458)^2*b^2*(1-b^2)^-1)/(1.60217657e-16); % kinematically limited max energy, in eV
k=6e-11*1e6; % g cm^-2 keV^-alpha_w -> kg mm^-2 keV^-alpha_w -> kg m^-2 keV^-alpha_w
density=1e3; % kg m^-3

Z=1;
Z_star=Z*(1-exp(-125*b*Z^(-2/3))); % effective charge of ion

energydensitylist=zeros(5,length(rlist));
for i=1:5
    for j=1:length(rlist)
        r=rlist(j);
        lowerlimit=(r*density/k)^(1/1.079); %keV
        if lowerlimit>1 %keV
            lowerlimit=(r.*density/k)^(1/1.667);
        end

        lowerlimit=lowerlimit*1000; %eV

        if lowerlimit<=capW % only run if range of energies required is less than max
            % perform integrals in keV units
            ruddintegral=@(W) ruddcs_integral(W,r,i,E_ion);
            wlist=lowerlimit:(capW-lowerlimit)/100000:capW;
            energydensitylist(i,j)=trapz(wlist,ruddintegral(wlist));
            energydensitylist(i,j)=Z_star.^2.*(1./(2.*pi.*r)).*energydensitylist(i,j); % eV/kg
            energydensitylist(i,j)=energydensitylist(i,j)*1.602e-19; % J/kg
        end
    end
end

energydensity=sum(energydensitylist(1:5,:),1);
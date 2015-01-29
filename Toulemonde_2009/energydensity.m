% Adapted from Waligorski, Hamm, and Katz (1986)
% http://digitalcommons.unl.edu/cgi/viewcontent.cgi?article=1099&context=physicskatz
% Calculates dose deposited in a coaxial cylindrical shell of thickness dr
% at a distance r from the path of an ion
% Last edited by Ping Lin/Gabriel, 27/01/15

function D2_r = energydensity(r,w,E_ion,Z)  
    % electron energy w in units of keV
    % ion energy E_ion in units of MeV
    switch nargin
        case 3
            Z=1; % charge of proton
    end
    c1=1.352817016; % Ne^4/mc^2 in units of keV mm^-1
    lorentzfactor=(E_ion/938)+1;
    b=(1-(1/lorentzfactor)^2)^0.5;
    if (w<1)
        alpha_w=1.079;
    else 
        alpha_w=1.667;
    end
    if (b<0.03)
        alpha_b=1.079;
    else
        alpha_b=1.667;
    end
    k=6e-6; % g cm^-2 keV^-alpha_w
    theta=range(k,w,alpha_w); % range of electron of energy w
    capW=2*9.1e-31*(3e8)^2*b^2*(1-b^2)^-1; % kinematically limited max energy
    T=range(k,capW,alpha_b); % max range
    Z_star=Z*(1-exp(-125*b*Z^(-2/3))); % effective charge of ion
    D1_r=zeros(1,length(r));
    D2_r=zeros(1,length(r));
    for i=1:length(r)
        D1_r(i) = c1*(Z_star^2)*(alpha_w*b^2*r(i)*(r(i)+theta))^-1*...
            (1-(r(i)+theta)/(r(i)+T))^(1/alpha_w);
        if (r(i)>0.1e-6) % units in mm NOT nm
            C=1.5e-6+5e-6*b;
            if (b<0.03)
                A=8*b^(1/3);
            else
                A=19*b^(1/3);
            end
            K_r=A*((r(i)-0.1e-6)/C)*exp(-((r(i)-0.1e-6)/C));
        else
            K_r=0;
        end
        D2_r(i) = D1_r(i)*(1+K_r);
    end
end
    
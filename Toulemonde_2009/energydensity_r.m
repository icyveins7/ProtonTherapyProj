% Adapted from Waligorski, Hamm, and Katz (1986)
% http://digitalcommons.unl.edu/cgi/viewcontent.cgi?article=1099&context=physicskatz
% Calculates dose deposited in a coaxial cylindrical shell of thickness dr
% at a distance r from the path of an ion
% Last edited by Ping Lin/Gabriel, 27/01/15
% 
% [D2_r,D1_r] = energydensity_r(r,E_ion)  
% [D2_r,D1_r] = energydensity_r(r,E_ion,w)  
% [D2_r,D1_r] = energydensity_r(r,E_ion,w,Z)

function [D2_r,D1_r] = energydensity_r(r,E_ion,varargin)  
    % electron energy w in units of keV
    % ion energy E_ion in units of MeV
    % r in units of mm
    switch nargin
        case 2 % default values
            w=0.010; %keV
            Z=1;
        case 3   
            w=varargin{1};
            Z=1; % charge of proton
        case 4
            w=varargin{1};
            Z=varargin{2};
        otherwise
            error('Invalid no. of args.');
    end
    density=1e-6;
    r=r*density; % mm -> kg/mm^2 by multiplying by density (kg/mm^3)
    c1=1.352817016; % Ne^4/mc^2 in units of keV mm^-1
    c1=c1*(1.60217657e-16); % J mm^-1
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
    k=6e-11; % g cm^-2 keV^-alpha_w -> kg mm^-2 keV^-alpha_w
    theta=range(k,w,alpha_w); % range of electron of energy w, kg mm^-2
    capW=(2*9.10938291e-31*(299792458)^2*b^2*(1-b^2)^-1)/(1.60217657e-16); % kinematically limited max energy, in keV
    T=range(k,capW,alpha_b); % max range, kg mm^-2
    Z_star=Z*(1-exp(-125*b*Z^(-2/3))); % effective charge of ion
    D1_r=zeros(1,length(r));
    D2_r=zeros(1,length(r));
    for i=1:length(r)
        if r(i)<=k*1^1.667
            alpha=1.079;
        else
            alpha=1.667;
        end
%         D1_r(i) = (c1*(Z_star^2)*(alpha*b^2*r(i)*(r(i)+theta))^(-1)*...
%             (1-(r(i)+theta)/(theta+T))^(1/alpha))*density; %J/kg, after multiplying by density (kg/mm^3), corrected?
        D1_r(i) = (c1*(Z_star^2)*(alpha_w*b^2*r(i)*(r(i)+theta))^(-1)*...
            (1-(r(i)+theta)/(theta+T))^(1/alpha_w))*density; %J/kg, after multiplying by density (kg/mm^3)
        if (r(i)>(0.1e-6*density)) % units in mm NOT nm, multiplied by the density to compare with modified r(i) value
            C=1.5e-6+5e-6*b;
            if (b<0.03)
                A=8*b^(1/3);
            else
                A=19*b^(1/3);
            end
            K_r=A*((r(i)/density-0.1e-6)/C)*exp(-((r(i)/density-0.1e-6)/C));
        else
            K_r=0;
        end
        % set all negative/complex numbers to 0
        if ~isreal(D1_r(i)) || real(D1_r(i))<0
            D1_r(i)=0;
        end       
        D2_r(i) = D1_r(i)*(1+K_r);
    end
end
    
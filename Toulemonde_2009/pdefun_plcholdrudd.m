function [c,f,s]=pdefun_plcholdrudd(x,t,u,dudx,rho,g,E_ion,bconst,S_e)
    % 1st equation electronic system
    % 2nd equation molecular system
    % rho in g/cm^3, E_ion in MeV, S_e in J/nm, x in nm, t in s, u in K
    
    % c=[3.0468e-25; 4.18*rho*1e-21]; %rho_e worked out by us, see next line
    % c=[3.053356543239642e-25; 4.18*rho*1e-21];% rho_e now KNOWN, calculated from equation in Waligorski paper
    specheatdensity=4.18;
    if u(2)>=373
        specheatdensity=2;
    end
    c=[1e-21; specheatdensity*rho*1e-21];
    f=[2e-7*dudx(1);6e-10*dudx(2)];
    if u(1)<=u(2)
        g=0;
    end

    E_r=rudd_density(x*1e-9,E_ion); % J kg^-1
    E_r=E_r*1e-3; %J g^-1
    E_t=energydensity_t(t);
    s=[-g*(u(1)-u(2)) + bconst.*S_e.*E_r.*E_t;
        g*(u(1)-u(2))];
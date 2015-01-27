% Adapted from Toulemonde, Surdutovich, and Solov'yov (2009)
% i-TS model: energy transfer into molecular subsystem
% Last edited by Ping Lin/Gabriel, 27/01/15

function dTdt = moleculartempfunc(T,t)
    rho=1; % density of liquid water, in g cm^-3
    C=4.18; % specific heat of water between 273-373K, in J g^-1 K^-1
    Ke=2; % electronic thermal conductivity, in W K^-1 cm^-1
    K=6e-3; % thermal conductivity, in W K^-1 cm^-1
    Te=310; % electron system temperature, in K
    lambda=2e-7; % mean free path, in cm
    g=Ke/(lambda^2); % electron-phonon coupling
    dTdt=(gradient(K*gradient(T))+g(Te-T))/(rho*C);
end


    
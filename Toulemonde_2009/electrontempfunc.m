% Adapted from Toulemonde, Surdutovich, and Solov'yov (2009)
% i-TS model: energy transfer into electron subsystem
% Last edited by Ping Lin/Gabriel, 27/01/15

function dTedt = electrontempfunc(Te,t,r,w,v,Z)
    Ce=1; % specific heat, in J K^-1 g^-1
    Ke=2; % thermal conductivity, in W K^-1 cm^-1
    T=310; % molecular system temperature, in K
    lambda=2e-7; % mean free path, in cm
    g=Ke/(lambda^2); % electron-phonon coupling
    dTedt=(gradient(Ke*gradient(Te))-g(Te-T)+energydensity(r,w,v,Z)*1e3)/Ce;
end
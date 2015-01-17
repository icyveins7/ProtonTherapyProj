% Equation 15 of Bortfeld paper

function D_z=dose(phi0,sigma,beta,alpha,gamma,E0,p,d,rho)
    steps=length(d);
    D_z=zeros(1,steps);
    R0=range(alpha,E0,p);
    z=zeta(R0,d,sigma);
    for i=1:steps
        D_z(i) = phi0*exp(-(z(i)^2)/4)*(sigma^(1/p))*gammafunction(1/p)*((sqrt(2*pi)*rho*p*alpha^(1/p)*(1+beta*R0))^-1)*((1/sigma)*yE_parabolic(0.5+1/p,-1*z(i))+(gamma*beta+beta/p)*yE_parabolic(0.5+1/p,-1*z(i)));
    end
end
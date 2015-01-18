function phi_z=fluence(phi0,beta,R0,z)
    phi_z=phi0*((1+beta*(R0-z))/(1+beta*R0));
end
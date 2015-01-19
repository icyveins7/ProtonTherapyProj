function D_nu=parafunc(nu,z)
    D_nu=2^(nu/2)*exp(-z.^2/4).*kummerU(-nu/2,1/2,0.5*z.^2);
end
function D_nu=parabolicfunc(nu,z)
    D_nu=2^(nu/2+1/4)*z.^(-1/2).*whittakerW(nu/2+1/4,-1/4,0.5*z.^2);
end
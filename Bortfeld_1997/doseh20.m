% Equation 29 of Bortfeld paper

function D_h20=doseh20(phi0,sigma,alpha,epsilon,E0,p,d)
    steps=length(d);
    D_h20=zeros(1,steps);
    R0=range(alpha,E0,p);
    z=zetafunc(R0,d,sigma);
    for i=1:steps
        D_h20(i) = phi0*exp(-((R0-z(i))^2)/(4*sigma^2))*sigma^0.565*((1+0.012*R0)^-1)*(11.26*sigma^-1*parabolicfunc(-0.565,-z(i))+(0.157+11.26*epsilon/R0)*parabolicfunc(-1.565,-z(i)));
    end
end
function x=xi(r,t,beta,rho1,epsilon)
    x=zeros(1,length(r));
    for i=1:length(r)
        x(i)= (r(i)/(beta*sqrt(t))) * (rho1/epsilon)^0.25;
    end
end
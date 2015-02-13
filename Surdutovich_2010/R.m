function Rad=R(beta,epsilon,rho1,t)
    Rad=zeros(1,length(t));
    for i=1:length(t)
        Rad(i)= beta * (t(i)^0.5) * (epsilon/rho1)^0.25;
    end
end
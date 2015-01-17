function z=zeta(R0,d,sigma)
    steps=length(d);
    z=zeros(1,steps);
    for i=1:steps
        z(i)=(R0-d)/sigma;
    end
end
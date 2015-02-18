function p=pressure_r(xi,G,Z,t_instant)
    p=zeros(1,length(xi));
    for i=1:length(xi)
        p(i)= 1.83 * (xi(i))^2 * G * Z / t_instant;
    end
end
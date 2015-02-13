function p2=pressure2(y,t)
    p2=zeros(1,length(t));
    for i=1:length(t)
        p2(i)= 1.83 * (2*y) * (y+1)^(-1) * t(i)^-1;
    end
end
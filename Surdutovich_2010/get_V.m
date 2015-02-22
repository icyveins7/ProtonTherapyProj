function vlist=get_V(y,xi)
    vlist=zeros(1,length(xi));
    for i=1:length(xi)
        syms v
        vlist(i) = vpasolve( (4/(y^2-1))*((y-1)/(y+1))^(1/y)*(((y*v-1)^(1-1/y))/(v*(2-y*v))) == xi(i)^2 , v, [0 100]);
    end
end
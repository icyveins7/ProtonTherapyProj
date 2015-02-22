function xi2=xi_square(y,V)
    xi2=zeros(1,length(V));
    for i=1:length(V)
        xi2(i)=(4/(y^2-1))*((y-1)/(y+1))^(1/y)*(((y*V(i)-1)^(1-1/y))/(V(i)*(2-y*V(i))));
    end
end
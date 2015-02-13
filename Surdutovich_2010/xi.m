function x=xi(r,R)
    x=zeros(1,length(r));
    for i=1:length(r)
        x(i)=r(i)/R(i);
    end
end
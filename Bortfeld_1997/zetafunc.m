% Last edited 18th Jan
% From Bortfeld/Schlegel, 1996

function z=zetafunc(R0,d,sigma)
    steps=length(d);
    z=zeros(1,steps);
    for i=1:steps
        z(i)=(R0-d(i))/sigma;
    end
end
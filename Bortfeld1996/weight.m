% Last edited: 16th Jan
% From Bortfeld/Schlegel, 1996

function W_R=weight(alpha,p,R,D0,rho,d_a,d_b)
    steps=length(R);
    W_R=zeros(1,steps);
    for i=1:steps
        if (R(i)<=d_b) && (R(i)>=d_a)
            W_R(i)=rho*D0*p*sin(pi/p)*alpha^(1/p)/(pi*(d_b-R(i))^(1/p));
        end
    end
end
% Last edited: 16th Jan
% From Bortfeld/Schlegel, 1996

function D_BP_d=depthdose(alpha,E0,p,d,rho)
    steps=length(d);
    D_BP_d=zeros(1,steps);
    R=range(alpha,E0,p);
    for i=1:steps
        if (d(i)<=R)
            D_BP_d(i) = (rho*p*alpha^(1/p)*(R-d(i))^(1-(1/p)))^-1;
        end
    end
end

% Last edited: 16th Jan
% From Bortfeld/Schlegel, 1996

function E_d=energy_dist(alpha,E0,p,d) 
    steps=length(d);
    E_d=zeros(1,steps);
    R=range(alpha,E0,p);
    for i=1:steps
        E_d(i)=((R-d(i))/alpha)^(1/p);
    end
end
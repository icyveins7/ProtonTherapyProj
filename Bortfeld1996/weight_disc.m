% Last edited: 16th Jan

function W_disc_R=weight_disc(alpha,p,D0,rho,d_a,d_b,itv)
    steps=floor(((d_b-d_a)/itv))+1;
    R=d_a:itv:d_b;
    W_disc_R=zeros(1,length(R));
    for i=1:steps
        if (R(i)==d_b)
            W_disc_R(i)=(rho*D0*p^2*alpha^(1/p)*sin(pi/p)/(pi*(p-1))) * (itv/2)^(1-1/p);
        else
            W_disc_R(i)=(rho*D0*p^2*alpha^(1/p)*sin(pi/p)/(pi*(p-1))) * ( (d_b-R(i)+itv/2)^(1-1/p) - (d_b - R(i) - itv/2)^(1-1/p) );
        end
    end
    W_disc_R=vertcat(W_disc_R,R);
end
    
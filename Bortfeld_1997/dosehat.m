% Last edited: 19th Jan
% Eq 11 of Bortfeld paper, dose with no straggling

function Dhat=dosehat(phi0,beta,alpha,gamma,E0,p,d,rho)
    Dhat=zeros(1,length(d));
    R0=range(alpha,E0,p);
    for i=1:length(d)
        if d(i)<=R0
            Dhat(i)=phi0*( ((R0-d(i))^(1/p-1)+(beta+gamma*beta*p)*(R0-d(i))^(1/p))/(rho*p*alpha^(1/p)*(1+beta*R0)) );
        else
            Dhat(i)=0;
        end
    end
end
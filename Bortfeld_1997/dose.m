% Equation 26 of Bortfeld paper
% toggle_approx=1 -> uses dosehat for z<R0-10sigma

function D_z=dose(phi0,sigma,beta,alpha,gamma,E0,p,d,rho,epsilon,toggle_approx)
    steps=length(d);
    D_z=zeros(1,steps);
    R0=range(alpha,E0,p);
    z=zetafunc(R0,d,sigma);
    for i=1:steps
        if (d(i)<=(R0-10*sigma)) && (toggle_approx==1)
            D_z(i)=dosehat(phi0,beta,alpha,gamma,E0,p,d(i),rho);
        else
            D_z(i) = phi0*exp(-(z(i)^2)/4)*(sigma^(1/p))*double(feval(symengine,'gamma',(1/p)))*...
                ((sqrt(2*pi)*rho*p*alpha^(1/p)*(1+beta*R0))^-1)*((1/sigma)*parafunc(-1/p,-z(i))+...
                (gamma*beta+beta/p+epsilon/R0)*parafunc(-1/p-1,-z(i)));
        end
    end
end
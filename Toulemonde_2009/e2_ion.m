function value=e2_ion(E)

% k = 1 to 5
Ej=[11.95, 14.7, 16.6, 33.3, 540]; % eV
yj=[12.5, 16.1, 19.4, 95, 220]; % eV
trij=[1.16, 1.31, 0.55, 1, 1]; % eV
% fj=[3.4, 3.1, 2.394, 1.595, 3.11]/10;
fj=[2.25, 2.06, 1.61, 1.21, 1.79]/10;

value=zeros(1,length(E));

for i=1:length(E)
    for j=1:5
        %make list of omegas to integrate over
        omegalist=Ej(j)-trij(j):2*trij(j)/1000:Ej(j)+trij(j);
        Dlist=fj(j).*yj(j).*E(i)./((omegalist.^2 - E(i).^2).^2 + yj(j).^2.*E(i).^2);
        stepfunc=omegalist<=E(i);
        gaussfunc=exp(-(omegalist-Ej(j)).^2/(2.*trij(j).^2));
        integrand=Dlist.*stepfunc.*gaussfunc;

        value(i)=value(i)+trapz(omegalist,integrand);
    end
end

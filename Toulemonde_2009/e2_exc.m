function value=e2_exc(E)

% k = 1 to 5
Ek=[8.17, 10.13, 11.31, 12.91, 14.5]; % eV
yk=[1.62, 2.2, 2.1, 3.1, 3.9]; % eV
fk=[0.118, 0.23, 0.1675, 0.285, 0.28]/10;

Dstar=zeros(5,length(E));
for k=1:5
    Dstar(k,:)=2.*fk(k).*yk(k).^3.*E.^3./( ( Ek(k).^2 - E.^2).^2 + yk(k).^2.*E.^2 ).^2;
end
value=sum(Dstar,1);

nu=-1/1.77;
d=17; %cm
E0=150; %E0 in units of MeV
alpha=2.2e-3;
p=1.77;
R0=range(alpha,E0,p);
sigma=0.012*(R0^0.935);
zc=(R0-d)/sigma;

steps=70;
x=1:steps;
y=zeros(1,steps);

for i=1:steps
    y(i)=exp(x(i)*zc-(x(i)^2)/2)*x(i)^(-nu-1);
end

figure(1);
plot(x,y);
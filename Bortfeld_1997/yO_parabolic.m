function [y,yp]=yO_parabolic(a,z)
%{
This function calculates the EVEN Weber functions (parabolic cylinder functions) 
which are solutions of the Weber equation (parabolic cylinder equation):

d^2y/dz^2+(z^2/4-a)y=0

where
  -Infinity<z<Infinity 
   a=real parameter

INPUTS:
-Infinity<z<Infinity  REAL independent variable (Vector or Matrix). 
a = real parameter

OUTPUTS:
y =  EVEN Weber function
yo = First derivative of the EVEN Weber function

For more information about Weber functions:

    Miguel A. Bandres and B.M. Rodriguez-Lara, "Nondiffracting accelerating waves: Weber waves and parabolic momentum"
    New Journal of Physics, 15(013054) (2013) ( http://goo.gl/XfaVqq )

    Miguel A. Bandres, J. C. Gutiérrez-Vega, and S. Chávez-Cerda, "Parabolic nondiffracting optical wavefields"
    Optics Letters, 29(1), 44-46 (2004) ( http://goo.gl/KhmqQY )

AUTHOR: Miguel A. Bandres
URL:    www.mabandres.com

%}

%% Internal Parameters
a=-a;                                         % to agree with our numerical method
[largo,ancho]=size(z); z=z(:)';      % Input Matrix -> Input Vector
index_paridad=-(z<0)+(z>=0); 
zmax=max(abs(z));
h=0.1;                                       % control number of analytic continuations
nca=ceil(zmax/h);                       % Number of analytic continuations 
[z,index]=sort(abs(z));                 % sort the input vector


%% Analytic Continuation
yi=[0,1];        % Boundary Condition
y=[]; yp=[];    % Initialized variables

for j=1:nca
    xo=(j-1)*h;
    xf=j*h;
    zi=z(find(xo<=z&z<xf)); 
    zi=[zi,xf];
    [yin,ypin]=yc_parabolic(yi,xo,a,zi);
    pint=length(yin)-1;  
    y=[y,yin(1:pint)]; 
    yp=[yp,ypin(1:pint)];       
    xo=xo+h;
    yi=[yin(pint+1),ypin(pint+1)]; 
end

if zmax==(nca*h)&(nca~=0)
    num=sum((z==zmax));
    y=[y,yin(pint+1)*ones(1,num)];
    yp=[yp,ypin(pint+1)*ones(1,num)];
end    
    
if nca~=0
    yv(index)=y;
    ypv(index)=yp;
    yv=yv.*index_paridad;
else
    yv=0;
    ypv=1;
end

%% Ouput
y=reshape(yv,[largo,ancho]);     % reshape output to original format
yp=reshape(ypv,[largo,ancho]);
    
return;

%% Subfuction to Calculate Polynomial Expansion
function [y,yp]=yc_parabolic(yi,xo,a,r)
N=100; % Number of polynomial expansion coefficients 

% Calculate coefficients
ac=parabolic_coe(yi,xo,a,N); 
acn=ac./gamma(1:N); 
acp=ac(2:N)./gamma(1:N-1); 

% Evaluate polynomials
y=polyval(acn(N:-1:1),(r-xo));
yp=polyval(acp(N-1:-1:1),(r-xo));
return;

%% Subfuction to Calculate coefficients of Polynomial Expansion 
function [yc]=parabolic_coe(yi,xo,a,N)
% Calculate the first N coefficientes of the polynomial expansion

yc=zeros(1,N);
yi=[yi,-((xo^2)/4+a)*yi(1)]; 
yc(1:3)=yi(1:3);
yc(4)=-(4*a+xo^2)/4*yc(2)-xo/2*yc(1);

for n=4:N-1
     yc(n+1)=-(4*a+xo^2)/4*yc(n+1-2)-(n-2)/2*xo*yc(n+1-3)-(n-2)*(n-3)/4*yc(n+1-4);
 end;

return;
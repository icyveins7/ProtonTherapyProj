function f = fcoeffunction(p,t,u,time,E_ion,bconst,S_e)

N = 2; % Number of equations
% Triangle point indices
it1 = t(1,:);
it2 = t(2,:);
it3 = t(3,:);

% Find centroids of triangles
xpts = (p(1,it1)+p(1,it2)+p(1,it3))/3;
ypts = (p(2,it1)+p(2,it2)+p(2,it3))/3;
rpts = (xpts.^2 + ypts.^2).^0.5;

% % dont think we need these
% [ux,uy] = pdegrad(p,t,u); % Approximate derivatives
% uintrp = pdeintrp(p,t,u); % Interpolated values at centroids 

nt = size(t,2); % Number of columns
f = zeros(N,nt); % Allocate f

% Now the particular functional form of f
% note that energydensity_r returns in J/kg, so must convert to J/g
% S_e should be in J/nm, input rpts as mm, E_ion in MeV
f(1,:) = bconst.*S_e.*energydensity_r(rpts.*1e-6,E_ion)./1000.*energydensity_t(time); % 1e25 factor???
f(2,:) = 0;
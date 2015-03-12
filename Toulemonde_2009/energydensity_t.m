% Adapted from Toulemonde (2006)
% Time dependence of dose density function
% Last edited by Ping Lin/Gabriel, 27/01/15
% 
% energydensity_t(t)
% energydensity_t(t,t0,s)

function D_t= energydensity_t(t,varargin) 
    switch nargin
        case 1 % default values
            t0=1e-15; % in seconds
            s=1e-14; % in seconds
        case 3
            t0=varargin{1};
            s=varargin{2};
        otherwise
            error('Invalid no. of args.');
    end
    
    D_t=zeros(1,length(t));
    
    for i=1:length(t)
        D_t(i)=exp(-(t(i)-t0)^2/(2*s^2));
    end
end
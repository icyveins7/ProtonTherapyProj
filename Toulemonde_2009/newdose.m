% testing angular distribution integration over proton path
function value=newdose(z1,r,z2,E0)
% z1,z2,r in m
% E0 in MeV

    % calculate angle phi from z1
    phi=atan((z2-z1)./r);
    theta=pi/2 - phi;
%     figure; plot(z1,phi); xlabel('z1'); ylabel('phi'); % tested
    
    % calc distances electrons from each z1 travel
    dists=sqrt((z1-z2).^2 + r.^2); %m
    density=1e3; % kg m^-3
    dists=dists.*density; % kg m^-2

    
    % Declaring constants/parameters
    alpha=2.2e-3;
    p=1.77;
    rho=1; %mass density of medium, g/cm^3
    R0=range(alpha,E0,p); %R in units of cm
    beta=0.012;
    gamma=0.6; %fraction of energy released in the nonelastic nuclear 
    %interactions that is absorbed locally
    phi0=1000; %primary particle fluence
    epsilon=0.1; %fraction of peak fluence in tail fluence
    LET_integral=@(z1) 100.*rho.*dose_C(phi0,beta,alpha,gamma,E0,p,z1.*100,rho,0,0,1)./fluence(phi0,beta,R0,z1.*100); % bortfeld takes in values in cm
        
    E_rem=zeros(1,length(z1));
    b=zeros(1,length(z1)); %list of beta values
    
%     to speed up evaluations we integrate only the first E_rem value with 'integral'
%     and following z1 values with trapz
    E_rem(1)=E0-integral(LET_integral,0,z1(1));
    LET=100.*rho.*dose_C(phi0,beta,alpha,gamma,E0,p,z1.*100,rho,0,0,1)./fluence(phi0,beta,R0,z1.*100);
    b(1)=sqrt( 1 - (938.272046/(E_rem(1)+938.272046 ))^2 );
    for i=2:length(E_rem)
%         E_rem(i)=E0-integral(LET_integral,0,z1(i)); % MeV
        E_rem(i)=E_rem(1)-trapz(z1(1:i),LET(1:i)); % MeV
        b(i)=sqrt( 1 - (938.272046/(E_rem(i)+938.272046 ))^2 );
    end
    
%     figure; plot(z1,E_rem); xlabel('z1'); ylabel('E_rem'); % tested
%     figure; plot(z1,b); xlabel('z1'); ylabel('beta'); % tested
    
    capW=(2.*9.10938291e-31.*(299792458).^2.*b.^2.*(1-b.^2).^-1)./(1.60217657e-19); % kinematically limited max energy, in eV
    w = capW.*(cos(theta)).^2; %eV
    
%     figure; plot(z1,theta); ylabel('theta'); % tested
%     figure; plot(z1,w); hold on; plot(z1,capW); legend('w','capW'); 
    
    n=3.343*10^28; % Number density of molecules for liquid water (molecules/m^3)
    Ne=n*10; % number density of electrons (10 per molecule) for liquid water (electrons/m^3)
    
    no_e=zeros(1,length(w));
    for j=1:length(w) % m^2
        no_e(j)=rudd_cs(w(j),1,E_rem(j))+rudd_cs(w(j),2,E_rem(j))+rudd_cs(w(j),3,E_rem(j))+rudd_cs(w(j),4,E_rem(j))+rudd_cs(w(j),5,E_rem(j));
    end
    no_e=0.5.*Ne.*no_e; % half factor
    
    difflist=diff(w);
    dw_dzlist=zeros(1,length(w));
    % for first and last use forward/backward approx
    dw_dzlist(1)=(w(2)-w(1))/(z1(2)-z1(1));
    dw_dzlist(end)=(w(end)-w(end-1))/(z1(end)-z1(end-1));
    for i=2:length(w)-1 % for all other values use central diff
        dw_dzlist(i)= (difflist(i)+difflist(i-1))/(z1(i+1)-z1(i-1));
    end
    
%     figure; plot(z1,dw_dzlist); ylabel('dw/dz');
    
    alpha_w=zeros(1,length(w));
    for l=1:length(w)
        if (w(l)<1000)
            alpha_w(l)=1.079;
        else 
            alpha_w(l)=1.667;
        end
    end
    k=6e-11*1e6; % g cm^-2 keV^-alpha_w -> kg mm^-2 keV^-alpha_w -> kg m^-2 keV^-alpha_w
    
    w_range=k.*(w./1000).^alpha_w;
    energy_per_e=zeros(1,length(w)); value = zeros(1,length(w));
    
%     figure; plot(z1,dists); hold on; plot(z1,w_range); 
%     legend('dists','w_range'); 
    
    for a=1:length(w_range)
        if w_range(a)>dists(a) % if range of electron is greater than distance from z1 to z2
            energy_per_e(a)=(w(a)./(alpha_w(a).*k.*(w(a)./1000).^(alpha_w(a)))) .* (1 - dists(a)./(k.*(w(a)./1000).^alpha_w(a))).^(1./alpha_w(a) - 1);
    
            value(a) = (2*pi.*r).^(-1) .* no_e(a) .* energy_per_e(a) .*abs(dw_dzlist(a));
            value(a) = value(a)*1.602e-19; % convert from eV/kg to J/kg
        else
            energy_per_e(a)=0;
            value(a)=0;
        end
    end
    
%     figure; plot(z1,energy_per_e);
%     figure; plot(z1,no_e);
    
    

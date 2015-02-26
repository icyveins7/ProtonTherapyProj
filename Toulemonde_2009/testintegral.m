function value=testintegral(w)
    % testing
    % w in keV
    density=1e-6;
    r=1e-6*density;
    
    I=0.078; %keV
%     I=0;

    if (w<1)
        alpha_w=1.079;
    else 
        alpha_w=1.667;
    end
    k=6e-11; % g cm^-2 keV^-alpha_w -> kg mm^-2 keV^-alpha_w
    
    value = 1./(k.*w.^(alpha_w-1)) .* (1 - r.*density./(k.*w.^alpha_w)).^(1./alpha_w - 1) .* (w+I).^(-2);
%     value = real(1./(k.*w.^(alpha_w-1)) .* (1 - r./(k.*w.^alpha_w)).^(1./alpha_w - 1) .* (w+I).^(-2));
%     value = imag(1./(k.*w.^(alpha_w-1)) .* (1 - r./(k.*w.^alpha_w)).^(1./alpha_w - 1) .* (w+I).^(-2));

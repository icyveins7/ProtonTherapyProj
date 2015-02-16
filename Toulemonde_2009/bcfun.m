function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t)
    pl=[0;0]; %ignored
    ql=[0;0]; %ignored
    pr=[ur(1)-310; ur(2)-310]; % at max radius, 310K boundary condition
    qr=[0;0];
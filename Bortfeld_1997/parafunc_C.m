% C library version of previous function, precision to 12 digits
% (equivalent to format long in Matlab
% Parabolic Cylinder Function
% Edited 21st Jan
% obtained from Table of Integrals, Series and Products by Gradshteyn,
% Ryzhik, 9.240, 2nd line

function [D_nu,terms]=parafunc_C(nu,z)
    D_nu=zeros(1,length(z));
    for i=1:length(z)
        terms=zeros(2,length(z));
        terms(1,i)=sqrt(pi)*calllib('Bortfeld_1997dll','F11',-nu/2,1/2,z(i)^2/2,12)/double(feval(symengine,'gamma',(1-nu)/2));
        terms(2,i)=sqrt(2*pi)*z(i)*calllib('Bortfeld_1997dll','F11',(1-nu)/2,3/2,z(i)^2/2,12)/double(feval(symengine,'gamma',-nu/2));
        D_nu(i)=2^(nu/2)*exp(-z(i)^2/4).*( terms(1,i) - terms(2,i) );
    end
end
function D_nu=parafunc2(nu,z)
%     term1=sqrt(pi)*double(feval(symengine,'hypergeom',-nu/2,1/2,z.^2/2))/double(feval(symengine,'gamma',(1-nu)/2));
%     term2=sqrt(2*pi)*z*double(feval(symengine,'hypergeom',(1-nu)/2,3/2,z.^2/2))/double(feval(symengine,'gamma',-nu/2));
    D_nu=zeros(1,length(z));
    for i=1:length(z)
        term1=sqrt(pi)*hypergeom(-nu/2,1/2,z(i)^2/2)/double(feval(symengine,'gamma',(1-nu)/2));
        term2=sqrt(2*pi)*z(i)*hypergeom((1-nu)/2,3/2,z(i)^2/2)/double(feval(symengine,'gamma',-nu/2));
        term1-term2
        D_nu(i)=2^(nu/2)*exp(-z(i)^2/4).*( term1 - term2 );
    end
end
% 9.246 Asymptotic Expansion for |z|>>1 and |z|>>|p|, first three terms
% Last edited 19th Jan
function [D_p,terms]=asymexp1(p,z)
    D_p=zeros(1,length(z));
    terms=zeros(2,length(z));
    for i=1:length(z)
        terms(1,i)=-p*(p-1)/(2*z(i)^2);
        terms(2,i)=p*(p-1)*(p-2)*(p-3)/(2*4*z(i)^4);
        D_p(i)=exp(-z(i)^2/4)*z(i)^p*(1+terms(1,i)+terms(2,i));
    end
end
    
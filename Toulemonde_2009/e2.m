function value = e2(E)
value=zeros(1,length(E));
for i=1:length(E)
%     value(i)=460.5316.*(e2_exc(E(i))+e2_ion(E(i)));
    value(i)=e2_C(E(i)); % C version
end
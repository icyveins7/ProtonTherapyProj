function value = e1 (E)

value=zeros(1,length(E));

for i=1:length(E)
    h=0.02;
    eprimelist=E(i):h:E(i)+25; % in front of singularity
    integrandlist=e2(eprimelist)./(eprimelist-E(i));
    value(i)=value(i)+trapz(eprimelist(2:end),integrandlist(2:end));
    eprimelist=E(i)-25:h:E(i); % behind singularity
    integrandlist=e2(eprimelist)./(eprimelist-E(i));
    value(i)=value(i)+trapz(eprimelist(1:end-1),integrandlist(1:end-1));
end

value=1+value./pi;

% % using C version, does not give smooth values
% for i=1:length(E)
%     value(i)=e1_C(E(i));
%     disp(i);
% end

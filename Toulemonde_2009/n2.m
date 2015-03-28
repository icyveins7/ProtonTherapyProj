function value = n2(E)

% call and save first to prevent double calls to functions
e1values=e1(E);
e2values=e2(E);

value = e2values./(e1values.^2 + e2values.^2);
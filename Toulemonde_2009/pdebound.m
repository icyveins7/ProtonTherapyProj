%  not needed?

function [qmatrix,gmatrix,hmatrix,rmatrix] = pdebound(p,e,u,time)

N = 2; % Set N = the number of equations
ne = size(e,2); % number of edges
qmatrix = zeros(N^2,ne);
gmatrix = zeros(N,ne);
hmatrix = zeros(N^2,2*ne);
rmatrix = zeros(N,2*ne);

for k = 1:ne
    x1 = p(1,e(1,k)); % x at first point in segment
    x2 = p(1,e(2,k)); % x at second point in segment
    xm = (x1 + x2)/2; % x at segment midpoint
    y1 = p(2,e(1,k)); % y at first point in segment
    y2 = p(2,e(2,k)); % y at second point in segment
    ym = (y1 + y2)/2; % y at segment midpoint
    switch e(5,k)
        case {some_edge_labels}
            % Fill in hmatrix,rmatrix or qmatrix,gmatrix
        case {another_list_of_edge_labels}
            % Fill in hmatrix,rmatrix or qmatrix,gmatrix
        otherwise
            % Fill in hmatrix,rmatrix or qmatrix,gmatrix
        
    end
end
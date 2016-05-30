function u = flowline(L,J,gamma,W,alpha,beta,V0)
% FLOWLINE  compute finite difference solution to elliptic PDE problem (i.e.
%   2-point BVP)
%     (W(x) u_x)_x - alpha(x) u = beta(x)     on  0 < x < L
%     u(0) = V0,   u_x(L) = gamma
% form:
%   u = flowline(L,J,gamma,W,alpha,beta,V0)
% where:
%   L     = length of domain
%   J     = number of subintervals
%   gamma = constant right-hand Neumann boundary condition
%   W     = vector of length J+1 (staggered grid values W_{j+1/2}, j=1:J+1)
%           (yes, last value should be at L+dx/2, past the end)
%   alpha = vector of length J+1 (regular grid values alpha_j, j=1:J+1)
%   beta  = same size as alpha
%   V0    = constant left-hand Dirichlet boundary condition
% returns:
%   u     = solution vector (regular grid values u_j, j=1:J+1)
% example: testflowline.m
% example: convanalysis.m
% main example: ssaflowline.m

dx = L / J;
rhs = dx^2 * beta(:); % a (J+1) length column vector
rhs(1) = V0;
rhs(J+1) = rhs(J+1) - 2 * gamma * dx * W(J+1);

A = sparse(J+1,J+1);  % allocates no space yet
A(1,1) = 1.0;
for j=2:J   % fill by rows
  A(j,j-1:j+1) = [ W(j-1), -(W(j-1) + W(j) + alpha(j) * dx^2), W(j) ];
end
A(J+1,J) = W(J) + W(J+1); 
A(J+1,J+1) = - (W(J) + W(J+1) + alpha(J+1) * dx^2);
%spy(A)   % uncomment to see that shape of A is tridiagonal (in small J cases)
%full(A)  % uncomment to see values (in small J cases)

% scale A by rows; otherwise Matlab says "close to singular or badly scaled."
scale = full(max(abs(A),[],2));    % column vector of row maximums
for j=1:J+1,  A(j,:) = A(j,:) ./ scale(j);  end
rhs = rhs ./ scale;

% solve by Matlab/Octave default methods
u = A \ rhs;

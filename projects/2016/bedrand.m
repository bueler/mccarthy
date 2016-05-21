function lamlims = bedrand(l,Nx,Ny,P,periodic)
% BEDRAND  Generate and plot P samples of 2d random bed topography on unit
% square [0,1]x[0,1], with Nx by Ny grid.  By default uses periodized
% squared-exponential covariance function
%     K(x,x') = exp(- (sin(pi(x1-x1')^2+sin(pi(x2-x2')^2) / (2*l^2))
% Otherwise (periodic=false), it is the classical one
%     K(x,x') = exp(- ((x1-x1')^2+(x2-x2')^2) / (2*l^2))
% with correlation length l.
% Examples:
%    >> bedrand(0.1,40,40)           # one sample, periodic
%    >> bedrand(1,40,40,4,false)     # four samples, smooth, not periodic
% Suggests convergence:
%    >> bedrand(.1,20,20)
%    >> bedrand(.1,40,40)
%    >> bedrand(.1,80,80)
% Probably o.k. for glacier use (with appropriate multiply and shift):
%    >> bedrand(.1,100,100,10)

if nargin < 4
    P = 1;
end
if nargin < 5
    periodic = true;
end

Lx = 1.0;
Ly = 1.0;
if periodic
    dx = Lx / Nx;
    dy = Ly / Ny;
    x = (0:dx:Lx-dx)';
    y = (0:dy:Ly-dy)';
else
    dx = Lx / (Nx-1);
    dy = Ly / (Ny-1);
    x = (0:dx:Lx)';
    y = (0:dy:Ly)';
end

[xx,yy] = ndgrid(x,y);
xxx = xx(:);
yyy = yy(:);

N = Nx * Ny;  % total number of points
Sigma = zeros(N,N);

for k = 1:Ny
    for j = 1:Nx
        if periodic
            d2 = sin(pi*abs(xxx-x(j))/Lx).^2 + sin(pi*abs(yyy-y(k))/Ly).^2;
        else
            d2 = (xxx - x(j)).^2 + (yyy - y(k)).^2;
        end
        Sigma(:,(k-1)*Nx + j) = exp(- d2 / (2.0*l*l) );
    end
end

Z = 400;
[V,D] = eigs(Sigma,Z);  % FIXME The parameter Z, which is the number of computed
                        % eigenvalues/vectors, is very important in both success
                        % (generating the right kind of variation) and speed.
                        % But I do not know how to set it properly.
Dhalf = sqrt(abs(diag(D)));

%figure(99)
%semilogy(1:Z,Dhalf), grid on

for m = 1:P
    Y = V * ( Dhalf .* (V' * randn(N,1)) );
    figure(m)
    imagesc(x,y,reshape(Y,Nx,Ny))
end


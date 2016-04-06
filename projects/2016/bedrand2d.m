function bedrand2d(l,Nx,Ny)
% BEDRAND2D  generate and plot one sample of 2d random bed topography on unit
% square [0,1]x[0,1], with Nx by Ny grid, using squared-exponential
% covariance function
%     K(x,x') = exp(- |x-x'|^2 / (2*l^2))
% with correlation length l where x = (x_1,x_2) and
%    |x-x'|^2 = (x_1-x_1')^2 + (x_2-x_2')^2
% examples:
%    >> bedrand2d(0.1,41,41)

Lx = 1.0;
Ly = 1.0;

N = Nx * Ny;  % total number of points

dx = Lx / (Nx-1);
dy = Ly / (Ny-1);
x = (0:dx:Lx)';
y = (0:dy:Ly)';
[xx,yy] = ndgrid(x,y);
xxx = xx(:);
yyy = yy(:);

Sigma = zeros(N,N);

for k = 1:Ny
  for j = 1:Nx
    d2 = (xxx - x(j)).^2 + (yyy - y(k)).^2;
    Sigma(:,(k-1)*Nx + j) = exp(- d2 / (2.0*l*l) );
  end
end

%d2 = (xxx - x(21)).^2 + (yyy - y(21)).^2;
%figure(3), surf(x,y,reshape(d2,Nx,Ny))

%sig = eig(Sigma);  
%figure(2),  semilogy(abs(sig))

%Shalf = sqrtm(Sigma);
%Y = Shalf * randn(N,1);

[V,D] = eigs(Sigma,40);
Dhalf = sqrt(diag(D));
Y = V * ( Dhalf .* (V' * randn(N,1)) );

figure(1),  surf(x,y,reshape(Y,Nx,Ny))


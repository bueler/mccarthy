function bedrand2d(l)
% BEDRAND2D  generate and plot one sample of 2d random bed topography on unit
% square [0,1]x[0,1] using squared-exponential covariance function
%     K(x,x') = exp(-|x-x'|^2 / (2*l^2))
% with correlation length l where x = (x_1,x_2) and
%    |x-x'|^2 = (x_1-x_1')^2 + (x_2-x_2')^2
% examples: >> bedrand2d(0.2)

Lx = 1.0;
Ly = 1.0;

Nx = 21;
Ny = 21;
N = Nx*Ny;

dx = Lx / (Nx-1);
dy = Ly / (Ny-1);
x = (0:dx:Lx)';
y = (0:dy:Ly)';

Sigma = zeros(N,N);

FIXME for k = 1:Ny
  for j = 1:Nx
    Sigma(:,(k-1)*Ny + j) = exp(- (x - x(j)).^2 / (2.0*l*l) );
end

%figure(2),  sig = eig(Sigma);  semilogy(abs(sig))

Y = sqrtm(Sigma) * randn(N,P);
plot(x,Y)


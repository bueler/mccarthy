function bedrand1d(l,P)
% BEDRAND1D  generate and plot P samples of 1d random bed topography on unit
% interval [0,1] using squared-exponential covariance function
%     K(x,x') = exp(-|x-x'|^2 / (2*l^2))
% with correlation length l
% examples: >> bedrand1d(0.2,5)   % smoother
%           >> bedrand1d(0.01,5)  % rougher

Lx = 1.0;  % length of bed
N = 201;   % number of points

dx = Lx / (N-1);
x = (0:dx:Lx)';

Sigma = zeros(N,N);

for j = 1:N
  Sigma(:,j) = exp(- (x - x(j)).^2 / (2.0*l*l) );
end

%figure(2),  sig = eig(Sigma);  semilogy(abs(sig))

Y = sqrtm(Sigma) * randn(N,P);
plot(x,Y)


function convanalysis(test)
% CONVANALYSIS  performs a convergence study of flowline.m using testflowline.m
% form:
%   convanalysis(test)
% where:
%   test = see comments in testflowline.m
% example: >> convanalysis
% notes:
%   1. execution time is a few seconds
%   2. use "more off" to flush Octave print statements
%   3. result for test=1 case meaningless; merely shows build-up of roundoff
%      error because solution is already exact for small J cases
%   3. results for test=2 and test=3 cases shows optimal convergence O(dx^2)

if nargin<1, test = 2; end

fprintf('convergence analysis of testflowline.m for test %d\n',test)
doublings = 8;
J = 10 * 2.^(0:doublings);   % =[10 20 40 ... 2560]
dx = 1.0 ./ J;               % domain has length one
for j=1:doublings+1
  err(j) = testflowline(test,J(j));
  fprintf('testflowline.m in J=%d (dx=%.3e) case gives error = %.4e\n',...
          J(j),dx(j),err(j));
end
pf = polyfit(log(dx),log(err),1);
loglog(dx,err,'*','markersize',12)
hold on, loglog(dx,exp(pf(1)*log(dx)+pf(2)),'r','linewidth',2), hold off
grid on,  xlabel('dx','fontsize',16),  ylabel('maximum error','fontsize',16)
result = sprintf('convergence rate O(dx^{%.5f})\n',pf(1));
disp(result),  text(2*dx(end-1),err(end-1),result,'fontsize',18)


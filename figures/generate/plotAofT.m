% PLOTAOFT 

gasConst_R = 8.31441;  % J/(mol K)
A_cold     = 3.61e-13; % Pa^-3 / s
A_warm     = 1.73e3;   % Pa^-3 / s
Q_cold     = 6.0e4;    % J / mol
Q_warm     = 13.9e4;   % J / mol
crit_temp  = 263.15;   % K
n = 3;

T = 223.15:0.1:273.15;
A = zeros(size(T));
A(T < crit_temp) = A_cold * exp( - Q_cold ./ (gasConst_R * T(T < crit_temp)) );
A(T >= crit_temp) = A_warm * exp( - Q_warm ./ (gasConst_R * T(T >= crit_temp)) );

max(A)/min(A)

semilogy(T-273.15,A,'linewidth',5.0)
xlabel('T (degrees C)','fontsize',24), ylabel('A(T)','fontsize',24)
text(-30,2e-26,'Paterson-Budd form for A(T)','fontsize',24)
grid on

% to create .eps in Octave:
%print -depsc2 ../pdffigs/AofT.eps


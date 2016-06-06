function plotSWdata
% Plot some of the data that went into
%     R. Sayag and M. G. Worster, 2013. "Axisymmetric gravity currents of
%     power-law fluids over a rigid horizontal surface," J. Fluid Mechanics
%     716, doi:10.1017/jfm.2012.545
% This is data for the constant-flux experiments.
% Code has additional Bueler edits to show all data.

%%% from email R. Sayag to E. Bueler, 5/20/13: %%%
% ... I attached the data (below) for the average front position that is
% plotted in figure 4b and a short .m file to plot it. Note that it includes
% the early-time front position (as in fig 4a) that is truncated from fig 4b
% since it is inconsistent with lubrication approximation in this paper. The
% pipe that supplies the fluid has inner diameter of 8-10 mm [LATER: was 8 mm]
% and I do agree that plug flow would be a very good approximation for the
% inside flow.

clf

d=10.3/547; %cm/px

Q = [3.8173 7.33 10.235];  % fluxes in gm/s; numbered (?) 232,250,270
files = {'rN_30rpm_V2B','rN_50rpm_V2','rN_70rpm_V2'};

style = {'kx','bo','r+'};
labels = {};
for k = 1:3
    load(files{k})
    t = (jRange-jRange(1)) * dt; % s
    rN = R * d; % cm
    hold on
    loglog(t(2:end),rN(2:end),style{k})   % t(1) = 0 so loglog removes anyway
    labels{k} = ['Q=',num2str(Q(k)),' gm/s'];
    xlabel('t [s]')
    ylabel('r_N [cm]')
end
hold off
legend(labels)


tic

PVol = [0:0.0001:0.01];
Vol = PVol*4*pi/3;
M_kill = zeros(size(Vol));
% M1 = M_kill;
% M2 = M1;

parfor i = 1:length(Vol)
[M_kill(i)] = EthanolDiff_BTCS(Vol(i),0.01,5);
% [M1(i)] = EthanolDiff_BTCS(Vol(i),0.2,5);
% [M2(i)] = EthanolDiff_CN(Vol(i),0.2,5);
end

plot(PVol,smooth(M_kill),'linewidth',2)
%plot(PVol,M_kill,PVol,M1,PVol,M2,'linewidth',2)
xlabel('Initial Volume')
ylabel('Percent Cells Killed')
alw = 0.75;    % AxesLineWidth
fsz = 14;      % Fontsize

toc
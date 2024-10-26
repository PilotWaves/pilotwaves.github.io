tic

PVol = 0.5;
Vol = PVol*4*pi/3;
U_thresh = [0.1:0.01:0.9];
t = [0:72/(length(U_thresh)-1):72];
U_t = exp(log(1.8)*t/72) - 0.9;

M_kill = zeros(size(U_thresh));
M_kill_t = M_kill;

parfor i = 1:length(U_thresh)
[M_kill(i)] = EthanolDiff_BTCS(Vol,U_thresh(i),5);
[M_kill_t(i)] = EthanolDiff_BTCS(Vol,U_t(i),5);
end

plot(U_thresh,smooth(M_kill),'linewidth',2)
hold on
plot(U_t,smooth(M_kill_t),'r--','linewidth',2)
hold off
xlabel('Threshold Concentration')
ylabel('Percent Cells Killed')
alw = 0.75;    % AxesLineWidth
fsz = 14;      % Fontsize
toc
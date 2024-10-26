tic

PVol = [0:.01:1];
Vol = PVol*4*pi/3;
U_thresh = 0.25;
eps = [1:0.1:10];

M_kill = zeros(size(Vol));

for j = 1:length(eps)
parfor i = 1:length(Vol)
[U,M_kill(i)] = EthanolDiff_Gen(Vol(i),U_thresh,eps(j));
end
plot(PVol,M_kill,'linewidth',2)
xlabel('Initial Average Concentration (Fraction of Tumor Volume)')
ylabel('Percent Cells Killed')
alw = 0.75;    % AxesLineWidth
fsz = 14;      % Fontsize
hold on
pause(0.1)
end

hold off

toc
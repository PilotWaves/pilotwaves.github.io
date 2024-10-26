tic

PVol = [0:.01:1];
Vol = PVol*4*pi/3;

a = [0:0.01:exp(1)-1];
a = fliplr(a);
U_thresh = 1 - log(a+1);

M_kill = zeros(length(Vol),length(U_thresh));

for j = 1:length(U_thresh)
parfor i = 1:length(Vol)
[M_kill(i,j)] = EthanolDiff_BTCS(Vol(i),U_thresh(j),5);
end
end

colormap(jet)
col = M_kill/100;
mesh(U_thresh,PVol,smoothdata(M_kill),col)
ylabel('Initial Average Concentration (Fraction of Tumor Volume)')
xlabel('Threshold (u_T)')
zlabel('Percent Cells Killed')
alw = 0.75;    % AxesLineWidth
fsz = 14;      % Fontsize

toc
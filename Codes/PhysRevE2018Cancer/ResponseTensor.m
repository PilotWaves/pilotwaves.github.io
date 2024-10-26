function [M_kill] = ResponseTensor

PVol = [0:.01:1];
Vol = PVol*4*pi/3;

U_thresh = [0:0.005:0.4995 0.5:0.01:1];

M_kill = zeros(length(Vol),length(U_thresh));

for j = 1:length(U_thresh)
parfor i = 1:length(Vol)
[M_kill(i,j)] = EthanolDiff_BTCS(Vol(i),U_thresh(j),5);
end
end


tic

PVol = [0:0.01:1];
Vol = PVol*4*pi/3;
U_thresh = [0:0.01:1];
LeastEffect = 0*U_thresh;

parfor j = 1:length(U_thresh)
for i = 1:length(Vol)
LE = LeastEffect_BTCS(Vol(i),U_thresh(j),5);
if LE ~= 0
    LeastEffect(j) = LE*3/(4*pi);
    break
end
end
end

plot(U_thresh,smooth(LeastEffect),'LineWidth',4)
xlabel('U_{thresh}')
ylabel('Least Effective Dose (Percent Tumor Size)')
alw = 0.75;    % AxesLineWidth
fsz = 14;      % Fontsize

toc
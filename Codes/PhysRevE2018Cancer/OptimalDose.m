tic

eps = [0.1:0.1:10];
OD = eps*0;
%slope = eps*0;
ut2 = 0.25;

parfor j = 1:length(eps)
    [OD2,~] = OptimalDose_Threshold(ut2,eps(j));
    OD(j) = OD2/(4*pi/3);
    %slope(j) = OD2/ut2/(4*pi/3);
end

% plot([0 U_thresh],[0 OD],'LineWidth',4)
% xlabel('U_{thresh} = a - b*exp(-c*t)')
% ylabel('Optimal Dose')
% alw = 0.75;    % AxesLineWidth
% fsz = 14;      % Fontsize

plot(eps,OD)

toc
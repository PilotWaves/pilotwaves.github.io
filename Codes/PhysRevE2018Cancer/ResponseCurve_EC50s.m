tic

clear all

%U_thresh = [1.9:0.01:2];

DataFull = [0.0572 0.0742 0.077; 0.0792 0.153 0.1411; 0.1227 0.2373 0.2676; 0.1298 0.3405 0.4242];
Data = DataFull(end-2:end,:);
DataVol = zeros(4,1);

t = [24 48 72];

MinDiff1 = 100*ones(3,1);
MinDiff2 = 100*ones(3,1);

for U_thresh = [0:0.01:3]
    MinDiff2(1) = abs(Data(end - 1,1) - EthanolDiff_BTCS(pi/6,U_thresh,5));
    MinDiff2(2) = abs(Data(end - 1,2) - EthanolDiff_BTCS(pi/6,U_thresh,5));
    MinDiff2(3) = abs(Data(end - 1,3) - EthanolDiff_BTCS(pi/6,U_thresh,5));
    
    if MinDiff2(1) <= MinDiff1(1)
       u1 = U_thresh;
       MinDiff1(1) = MinDiff2(1);
    else
        break
    end
    if MinDiff2(2) <= MinDiff1(2)
       u2 = U_thresh;
       MinDiff1(2) = MinDiff2(2); 
    end
    if MinDiff2(3) <= MinDiff1(3)
       u3 = U_thresh;
       MinDiff1(3) = MinDiff2(3); 
    end
end

u0 = (u1*u3 - u2^2)/(u1 + u3 - 2*u2);
c = log((u0-u1)/(u0-u2))/24;
b = (u0 - u1)^2/(u0-u2);

U_thresh = u0 - b*exp(-c*t);

PVol = [0:.01:0.25];
Vol = PVol*4*pi/3;

M_kill = zeros(length(Vol),length(U_thresh));

for j = 1:length(U_thresh)
parfor i = 1:length(Vol)
[M_kill(i,j)] = EthanolDiff_BTCS(Vol(i),U_thresh(j),5);
end
M_kill(:,j) = smooth(M_kill(:,j));
end

AvDiff = ones(length(Vol),length(DataFull(:,3)))*10^6;

for i = 1:length(PVol)
    for j = 1:length(DataFull(:,3))
       AvDiff(i,j) = abs(DataFull(j,3) - M_kill(i,3));
       if AvDiff(i,j) == min(AvDiff(:,j))
        DataVol(j) = PVol(i);
       end
    end
end

set(0,'defaultAxesFontSize',14)

ec50_bootstrap([0; DataVol]',[0; DataFull(:,3)]');
hold on

plot(PVol,M_kill(:,3),'k','linewidth',4)
xlabel('Initial Average Concentration (Fraction of Tumor Volume)')
ylabel('Percent Cells Killed')
alw = 0.75;    % AxesLineWidth
axis([0 max(DataVol) 0 max(DataFull(end,:))])

%  txt = ['[',num2str(quant(1)),',',num2str(quant(3)),']'];
%  annotation('textbox',[0.55 0.15 0.1 0.1],'FontSize',14,'String',txt,'FitBoxToText','on');

% DataVolSpline = [DataVol(1):0.1:DataVol(3)];
% Data1 = spline(DataVol,Data(:,1),DataVolSpline);
% Data2 = spline(DataVol,Data(:,2),DataVolSpline);
% Data3 = spline(DataVol,Data(:,3),DataVolSpline);
% 
% plot(PVol,M_kill(:,1),PVol,M_kill(:,2),PVol,M_kill(:,3),DataVolSpline,Data1,'k--',...
%     DataVolSpline,Data2,'r--',DataVolSpline,Data3,'b--','linewidth',2)

toc
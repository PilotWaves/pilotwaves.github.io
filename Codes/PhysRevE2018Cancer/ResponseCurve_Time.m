tic

clear all

%U_thresh = [1.9:0.01:2];

Data = [0.2246 0.1894 0.2745; 0.2404 0.3802 0.4735; 0.2708 0.5238 0.5417];
DataVol = zeros(3,1);

t = [24 48 72];

MinDiff1 = 100*ones(3,1);
MinDiff2 = 100*ones(3,1);

for U_thresh = [0:0.01:3]
    MinDiff2(1) = abs(Data(2,1) - EthanolDiff_BTCS(pi/6,U_thresh,5));
    MinDiff2(2) = abs(Data(2,2) - EthanolDiff_BTCS(pi/6,U_thresh,5));
    MinDiff2(3) = abs(Data(2,3) - EthanolDiff_BTCS(pi/6,U_thresh,5));
    
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

AvDiff = ones(size(M_kill))*10^6;

for i = 1:length(PVol)
    if (abs(Data(1,1) - M_kill(i,1)) + abs(Data(1,2) - M_kill(i,2)) < ...
       abs(Data(1,1) - M_kill(i,1)) + abs(Data(1,3) - M_kill(i,3)) &&...
       abs(Data(1,1) - M_kill(i,1)) + abs(Data(1,2) - M_kill(i,2)) < ...
       abs(Data(1,3) - M_kill(i,3)) + abs(Data(1,2) - M_kill(i,2)))
    AvDiff(i,1) = abs(Data(1,1) - M_kill(i,1)) + abs(Data(1,2) - M_kill(i,2));
    else if (abs(Data(1,1) - M_kill(i,1)) + abs(Data(1,3) - M_kill(i,3)) < ...
       abs(Data(1,1) - M_kill(i,1)) + abs(Data(1,2) - M_kill(i,2)) &&...
       abs(Data(1,1) - M_kill(i,1)) + abs(Data(1,3) - M_kill(i,3)) < ...
       abs(Data(1,3) - M_kill(i,3)) + abs(Data(1,2) - M_kill(i,2)))
    AvDiff(i,1) = abs(Data(1,1) - M_kill(i,1)) + abs(Data(1,3) - M_kill(i,3));
        else
    AvDiff(i,1) = abs(Data(1,3) - M_kill(i,3)) + abs(Data(1,2) - M_kill(i,2));
        end
    end
    
    if (abs(Data(2,1) - M_kill(i,1)) + abs(Data(2,2) - M_kill(i,2)) < ...
       abs(Data(2,1) - M_kill(i,1)) + abs(Data(2,3) - M_kill(i,3)) &&...
       abs(Data(2,1) - M_kill(i,1)) + abs(Data(2,2) - M_kill(i,2)) < ...
       abs(Data(2,3) - M_kill(i,3)) + abs(Data(2,2) - M_kill(i,2)))
    AvDiff(i,2) = abs(Data(2,1) - M_kill(i,1)) + abs(Data(2,2) - M_kill(i,2));
    else if (abs(Data(2,1) - M_kill(i,1)) + abs(Data(2,3) - M_kill(i,3)) < ...
       abs(Data(2,1) - M_kill(i,1)) + abs(Data(2,2) - M_kill(i,2)) &&...
       abs(Data(2,1) - M_kill(i,1)) + abs(Data(2,3) - M_kill(i,3)) < ...
       abs(Data(2,3) - M_kill(i,3)) + abs(Data(2,2) - M_kill(i,2)))
    AvDiff(i,2) = abs(Data(2,1) - M_kill(i,1)) + abs(Data(2,3) - M_kill(i,3));
        else
    AvDiff(i,2) = abs(Data(2,3) - M_kill(i,3)) + abs(Data(2,2) - M_kill(i,2));
        end
    end
    
    if (abs(Data(3,1) - M_kill(i,1)) + abs(Data(3,2) - M_kill(i,2)) < ...
       abs(Data(3,1) - M_kill(i,1)) + abs(Data(3,3) - M_kill(i,3)) &&...
       abs(Data(3,1) - M_kill(i,1)) + abs(Data(3,2) - M_kill(i,2)) < ...
       abs(Data(3,3) - M_kill(i,3)) + abs(Data(3,2) - M_kill(i,2)))
    AvDiff(i,3) = abs(Data(3,1) - M_kill(i,1)) + abs(Data(3,2) - M_kill(i,2));
    else if (abs(Data(3,1) - M_kill(i,1)) + abs(Data(3,3) - M_kill(i,3)) < ...
       abs(Data(3,1) - M_kill(i,1)) + abs(Data(3,2) - M_kill(i,2)) &&...
       abs(Data(3,1) - M_kill(i,1)) + abs(Data(3,3) - M_kill(i,3)) < ...
       abs(Data(3,3) - M_kill(i,3)) + abs(Data(3,2) - M_kill(i,2)))
    AvDiff(i,3) = abs(Data(3,1) - M_kill(i,1)) + abs(Data(3,3) - M_kill(i,3));
        else
    AvDiff(i,3) = abs(Data(3,3) - M_kill(i,3)) + abs(Data(3,2) - M_kill(i,2));
        end
    end
    
    if AvDiff (i,1) == min(AvDiff(:,1))
        DataVol(1) = PVol(i);
    end
    if AvDiff (i,2) == min(AvDiff(:,2))
        DataVol(2) = PVol(i);
    end
    if AvDiff (i,3) == min(AvDiff(:,3))
        DataVol(3) = PVol(i);
    end
end

set(0,'defaultAxesFontSize',14)

%LED = U_thresh(3)*1.1990/(4*pi/3);

%quant = ec50_bootstrap([LED; DataVol]',[0; Data(:,3)]');

plot(PVol,M_kill(:,1),PVol,M_kill(:,2),PVol,M_kill(:,3),'k',DataVol,Data(:,1),'bp',...
    DataVol,Data(:,2),'rp',DataVol,Data(:,3),'kp','linewidth',4)
% hold on
% errorbar(quant(2),0.5,quant(1)-quant(2),quant(3)-quant(2),'horizontal','gs','MarkerSize',15,'linewidth',2)
% hold off
xlabel('Initial Average Concentration (Fraction of Tumor Volume)')
ylabel('Percent Cells Killed')
alw = 0.75;    % AxesLineWidth

% DataVolSpline = [DataVol(1):0.1:DataVol(3)];
% Data1 = spline(DataVol,Data(:,1),DataVolSpline);
% Data2 = spline(DataVol,Data(:,2),DataVolSpline);
% Data3 = spline(DataVol,Data(:,3),DataVolSpline);
% 
% plot(PVol,M_kill(:,1),PVol,M_kill(:,2),PVol,M_kill(:,3),DataVolSpline,Data1,'k--',...
%     DataVolSpline,Data2,'r--',DataVolSpline,Data3,'b--','linewidth',2)

toc
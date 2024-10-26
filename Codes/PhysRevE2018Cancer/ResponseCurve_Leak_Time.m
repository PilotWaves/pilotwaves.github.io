tic

clear all

%U_thresh = [1.9:0.01:2];

Data = [0.0118 0.0203 0.0378; 0.0206 0.0385 0.1033; 0.0379 0.0682 0.193];
DataDose = [0.316, 1, 3.16];

dr = 0.01;
r = [0:dr:1];
width = 1;
r1 = round(length(r)*width);
height = (pi/6)/(sum((r(1:r1 - 1).^2).*exp(1 - 1./(1 - (r(1:r1 - 1)/r(r1)).^2)))*dr*4*pi);

MinDiff1 = 100*ones(3,1);
MinDiff2 = 100*ones(3,1);

for U_thresh = [0:0.001:height]
    MinDiff2(1) = abs(Data(1,3) - EthanolDiff_BTCS(pi/6,U_thresh,5));
    MinDiff2(2) = abs(Data(2,3) - EthanolDiff_BTCS(pi/6,U_thresh,5));
    MinDiff2(3) = abs(Data(3,3) - EthanolDiff_BTCS(pi/6,U_thresh,5));
    %MinDiff2(4) = abs(Data(4,3) - EthanolDiff_BTCS(pi/6,U_thresh,5));
    
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

%theta = log((u2 - height)*(u2 - u1)/((u3 - height)*(u1 - u3)))/log(DataDose(3)/DataDose(2));
%EC50 = (((u3 - height)*DataDose(3)^theta - (u2 - height)*DataDose(2)^theta)/(u2 - u3))^(1/theta);
%u_max = real(u3 - (height - u3)/(EC50/DataDose(3))^theta);
u_max = u3;

MinDiff_u_max = 100;
Resp_u3 = EthanolDiff_BTCS(pi/6,u3,5);
for i = 1:1000
    u_max = u_max - (i-1)/1000;
F1 = log((height - u1)/(u1 - u_max)); 
F2 = log((height - u2)/(u2 - u_max));  

EC50 = exp((F1*log(DataDose(2)) - F2*log(DataDose(1)))/(F1 - F2));
theta = F2/log(EC50/DataDose(2));

% EC50 = exp((log(height*u3 - 1)*log(DataDose(2)) - log(height/u2 - 1)*log(DataDose(3)))/log((height/u3 - 1)/(height/u2 - 1)));
% theta = log(height/u3 - 1)/log(EC50/DataDose(3));

U_thresh_72 = u_max + (height-u_max)./(1 + (EC50./DataDose(3)).^theta);
Resp_test = EthanolDiff_BTCS(pi/6,U_thresh_72,5);
if abs(Resp_u3 - Resp_test) <= MinDiff_u_max
    MinDiff_u_max = abs(Resp_u3 - Resp_test);
else
    u_max = u_max + (i-1)/1000;
    F1 = log((height - u1)/(u1 - u_max)); 
F2 = log((height - u2)/(u2 - u_max));  

EC50 = exp((F1*log(DataDose(2)) - F2*log(DataDose(1)))/(F1 - F2));
theta = F2/log(EC50/DataDose(2));
    break
end
end

d = exp(-4.6 + [0:0.1:5]*1.15);
U_thresh_72 = u_max + (height-u_max)./(1 + (EC50./d).^theta);

MinDiff1 = 100*ones(3,1);
MinDiff2 = 100*ones(3,1);

for U_thresh = [0:0.001:height]
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

U_thresh_48 = U_thresh_72*u2/u3;
U_thresh_24 = U_thresh_72*u1/u3;

M_kill_72 = 0*length(U_thresh_72);
M_kill_48 = M_kill_72;
M_kill_24 = M_kill_72;

parfor j = 1:length(U_thresh_72)
M_kill_72(j) = EthanolDiff_BTCS(pi/6,U_thresh_72(j),5);
M_kill_48(j) = EthanolDiff_BTCS(pi/6,U_thresh_48(j),5);
M_kill_24(j) = EthanolDiff_BTCS(pi/6,U_thresh_24(j),5);
end

M_kill_72 = smooth([0 M_kill_72]);
M_kill_48 = smooth([0 M_kill_48]);
M_kill_24 = smooth([0 M_kill_24]);

set(0,'defaultAxesFontSize',14)

%LED = U_thresh(3)*1.1990/(4*pi);

%quant = ec50_bootstrap([LED; DataVol]',[0; Data(:,3)]');

d = [0 d];
plot(d,M_kill_24,d,M_kill_48,d,M_kill_72,'k',DataDose,Data(:,1),'bp',...
    DataDose,Data(:,2),'rp',DataDose,Data(:,3),'kp','linewidth',4)
% hold on
% errorbar(quant(2),0.5,quant(1)-quant(2),quant(3)-quant(2),'horizontal','gs','MarkerSize',15,'linewidth',2)
% hold off
xlabel('Molarity (\mu M)')
ylabel('Percent Cells Killed')
alw = 0.75;    % AxesLineWidth
axis([0 3.16 0 1.1*max(Data(:,3))])

% DataVolSpline = [DataVol(1):0.1:DataVol(3)];
% Data1 = spline(DataVol,Data(:,1),DataVolSpline);
% Data2 = spline(DataVol,Data(:,2),DataVolSpline);
% Data3 = spline(DataVol,Data(:,3),DataVolSpline);
% 
% plot(leak,M_kill(:,1),leak,M_kill(:,2),leak,M_kill(:,3),DataVolSpline,Data1,'k--',...
%     DataVolSpline,Data2,'r--',DataVolSpline,Data3,'b--','linewidth',2)

toc
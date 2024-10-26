function [ED50,R50] = EC50(Data,R_mid)

t = [24 48 72];

MinDiff1 = 100*ones(3,1);
MinDiff2 = 100*ones(3,1);

for U_thresh = [0:0.001:3]
    MinDiff2(1) = abs(Data(1) - EthanolDiff_BTCS(pi/6,U_thresh,5));
    MinDiff2(2) = abs(Data(2) - EthanolDiff_BTCS(pi/6,U_thresh,5));
    MinDiff2(3) = abs(Data(3) - EthanolDiff_BTCS(pi/6,U_thresh,5));
    
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
MidVol = pi/3;

R50 = EthanolDiff_BTCS(MidVol,U_thresh(3),5);

if abs(R_mid - R50) < 0.001
       ED50 = MidVol;
        
else if R50 > R_mid
    left = 0; right = MidVol;
    mid = (left + right)/2;
for i = 1:100000
    R50 = EthanolDiff_BTCS(mid,U_thresh(3),5);
    
    if abs(R_mid - R50) < 0.001
        ED50 = mid;
       break
    end
    
    if R50 < R_mid
        left = mid;
    else
        right = mid;
    end
    mid = (left + right)/2;
end
ED50 = mid;
    else
            left = MidVol; right = 4*(4*pi/3);
    mid = (left + right)/2;
for i = 1:100000
    R50 = EthanolDiff_BTCS(mid,U_thresh(3),5);
    
    if abs(R_mid - R50) < 0.001
        ED50 = mid;
       break
    end
    
    if R50 < R_mid
        left = mid;
    else
        right = mid;
    end
    mid = (left + right)/2;
end
ED50 = mid;
    end
end
function [ED100,R100] = OptimalDose_Threshold(U_thresh,eps)

R100 = EthanolDiff_BTCS(10*U_thresh,U_thresh,eps);

if abs(0.99 - R100) < 1e-3
       ED100 = 10*U_thresh;
        
else if R100 > 0.99
    left = 0; right = 10*U_thresh;
    mid = (left + right)/2;
for i = 1:100000
    R100 = EthanolDiff_BTCS(mid,U_thresh,eps);
    
    if abs(0.99 - R100) < 1e-3
        ED100 = mid;
       break
    end
    
    if R100 < 0.99
        left = mid;
    else
        right = mid;
    end
    mid = (left + right)/2;
end
ED100 = mid;
    else
            left = 10*U_thresh; right = 30*U_thresh;
    mid = (left + right)/2;
for i = 1:100000
    R100 = EthanolDiff_BTCS(mid,U_thresh,eps);
    
    if abs(0.99 - R100) < 1e-3
        ED100 = mid;
       break
    end
    
    if R100 < 0.99
        left = mid;
    else
        right = mid;
    end
    mid = (left + right)/2;
end
ED100 = mid;
    end
end
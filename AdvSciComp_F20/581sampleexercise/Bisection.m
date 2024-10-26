function [A12, A13_2] = Bisection(xa, xb)

    xmid = (xa+xb)/2; A12 = xmid; %f(xa) > 0 and f(xb) < 0
    while  abs(xmid - xa) > 1e-1
        if -xmid - cos(xmid) > 0
            xa = xmid;
            else if -xmid - cos(xmid) < 0
            xb = xmid;
            else
                A12(length(A12)+1) = xmid;
                break
            end
        end
        xmid = (xa+xb)/2;
        A12(length(A12)+1) = xmid; %Solution to Problem 3b
    end
    A12 = A12';
    A13_2 = length(A12);
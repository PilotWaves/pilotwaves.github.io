function [A11, A13_1] = Newton_Ralphson(x)

    A11 = -3;
    while abs(-(-x-cos(x))/(-1+sin(x))) > 1e-6
        x =  x - (-x-cos(x))/(-1+sin(x)); % x_n+1 = x_n - f(x_n)/f'(x_n)
        A11(length(A11)+1) = x; %Solution to Problem 3a
    end
    A11 = A11(1:end)';
    A13_1 = length(A11)-1;
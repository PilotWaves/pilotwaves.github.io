function[out] = coupledFH(t,in_vec,param_vec)
v1 = in_vec(1);
w1 = in_vec(2);
v2 = in_vec(3);
w2 = in_vec(4);

a1 = param_vec(1);
a2 = param_vec(2);
b = param_vec(3);
c = param_vec(4);
I = param_vec(5);
d12 = param_vec(6);
d21 = param_vec(7);

out = [-v1^3 + (1+a1)*v1^2 - a1*v1 - w1 + I + d12*v2; ...
    b*v1 - c*w1; ...
    -v2^3 + (1 + a2)*v2^2 - a2*v2 - w2 + I + d21*v1; ...
    b*v2 - c*w2];
end
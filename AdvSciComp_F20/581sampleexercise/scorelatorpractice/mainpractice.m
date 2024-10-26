clear all; close all; clc

% exercise 1 - build a matrix

A=[34 45; 17 6];
save A1.dat A -ascii


% exercise 2 -  matrix operations

A = [1 2;-1 1];
B = [2 0;0 2];
C = [2 0 -3;0 0 -1];
D = [1 2;2 3;-1 0];
x = [1;0];
y = [0;1];
z = [1;2;-1];


A2 = A+B;
A3 = 3*x-4*y;
A4 = A*x;
A5 = B*(x-y);
A6 = D*x;
A7 = D*y+z;
A8 = A*B;
A9 = B*C;
A10 = C*D;

save A2.dat A2 -ascii
save A3.dat A3 -ascii
save A4.dat A4 -ascii
save A5.dat A5 -ascii
save A6.dat A6 -ascii
save A7.dat A7 -ascii
save A8.dat A8 -ascii
save A9.dat A9 -ascii
save A10.dat A10 -ascii

% exercise 3 -  rootfinding
[appvalsNR,niterNR]=newtonraphson(-3,10^-6);
[appvalsBis,niterBis]=bisection(-3,1,10^-6);

A11 = appvalsNR;
A12 = appvalsBis;
A13 = [niterNR niterBis];

save A11.dat A11 -ascii
save A12.dat A12 -ascii
save A13.dat A13 -ascii



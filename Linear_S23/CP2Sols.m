A = [1.1,.2,-.2,.5;
     .2,.9,.5,.3;
     .1,0.,1.,.4;
     .1,.1,.1,1.2];
b = [1;0;1;0];

L = tril(A,-1);
U = triu(A,+1);
D = diag(A);
LpD = tril(A); %L + D
LpU = L + U; % L + U
c_gs = Forsub([LpD,b]);
c_j = b./D;

y_gs = zeros(length(b),1);
y_j = zeros(length(b),1);

E_gs = zeros(4,1);
E_j = E_gs;
T_gs = E_gs;
T_j = E_gs;

count_j = 1;
count_gs = 1;

for i = 1:100000
    y_j = -(LpU*y_j)./D + c_j;
    y_gs = U*y_gs;
    y_gs = -Forsub([LpD,y_gs]) + c_gs;
    
   if count_j < 5 && max(abs(A*y_j - b)) < (1e-2)^count_j && max(abs(A*y_j - b)) > E_j(count_j)
       T_j(count_j) = i; E_j(count_j) = max(abs(A*y_j - b));
       count_j = count_j + 1;
   end
   
   if count_gs < 5 && max(abs(A*y_gs - b)) < (1e-2)^count_gs && max(abs(A*y_gs - b)) > E_gs(count_gs)
       T_gs(count_gs) = i; E_gs(count_gs) = max(abs(A*y_gs - b));
       count_gs = count_gs + 1;
   end
   
   if count_j > 4 && count_gs > 4
       break
   end
   
end

A1 = T_j(1); A2 = E_j(1); A3 = T_j(2); A4 = E_j(2); A5 = T_j(3); A6 = E_j(3); A7 = T_j(4); A8 = E_j(4); A9 = T_gs(1); A10 = E_gs(1); A11 = T_gs(2); A12 = E_gs(2); A13 = T_gs(3); A14 = E_gs(3); A15 = T_gs(4); A16 = E_gs(4);


x = [ .9;    % S
           .09;   % I
           .01 ]; % R
       
p = 0;
       
M = eye(3,3) + [-1/200 - p, 0, 1/10000; 1/200, -1/1000, 0; p, 1/1000, -1/10000];

D0 = 0;
F0 = 0;
for i = 1:100000
   x = M*x;
   if x(2) > 0.5 && D0 == 0
       D0 = i;
   end
   if abs(F0-x(2)) < 1e-8
      break 
   end
   F0 = x(2);
end

x = [ .9;    % S
           .09;   % I
           .01 ]; % R
p = 2/1000;
M = eye(3,3) + [-1/200 - p, 0, 1/10000; 1/200, -1/1000, 0; p, 1/1000, -1/10000];
D1 = 0;
F1 = 0;
for i = 1:100000
   x = M*x;
   if x(2) > 0.5 && D1 == 0
       D1 = i;
   end
   if abs(F1-x(2)) < 1e-8
      break 
   end
   F1 = x(2);
end

A17 = D0; A18 = F0; A19 = D1; A20 = F1;



M_k = zeros(22,21);
coeffs = [1:2:41];
coeffs = 1./coeffs;

M_k(1,:) = coeffs;
M_k(2:end,:) = diag(coeffs);
A21 = M_k;

a = 0*coeffs;
for i = 1:length(a)
   a(i) = (-1)^(i-1)/factorial(i-1); 
end

A22 = sum(M_k*a');



%save('CP2Sols.mat','A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','A13','A14','A15','A16','A17','A18','A19','A20','A21','A22');
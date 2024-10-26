%%%  Problem 1 %%%

A = [1, 0.1, 1, 0.1, 1, 1; 0.5, 0, 0.5, 0.4, 0.5, 1; 0.1, 0.1, 0.1, 0.1, 0.1, 0.1; 0, 1, 0, 0, 1, 1];
A1 = A;

b = [5; 2.5; 1.5; 5]; A2 = b; A3 = A*A';

x = A'*((A*A')\b); A4 = x;

for i = 1:length(x)
xifloor = x; xiceil = x; xifloor(i) = floor(xifloor(i)); xiceil(i) = ceil(xiceil(i));
if norm(A*xifloor - b) > norm(A*xiceil - b)
x(i) = ceil(x(i));
else
x(i) = floor(x(i));
end
end

A5 = x;

%%%  Problem 2 %%%

A6 = 17;

%%%  Problem 3 %%%

r = RandStream('mt19937ar','Seed',1234);
A = diag([1:1000]) + .01*r.randn(1000);
b = ones(length(A),1);

y= zeros(length(b),1);

E = 1e-8;
Tj = 0; Tgs = 0;

while max(abs(A*y - b)) > E
    Tj = Tj+1;
    y = ParallelJacobi(y,A,b);
end
Ej = max(abs(A*y - b));
A7 = y; A9 = Tj; A10 = Ej;

y= zeros(length(b),1);

while max(abs(A*y - b)) > E
    Tgs = Tgs+1;
    y = LinearGaussSeidel(y,A,b);
end
Egs = max(abs(A*y - b));

A8 = y; A11 = Tgs; A12 = Egs;

save('CP_Bonus_Sols.mat','A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12');
r = RandStream('mt19937ar','Seed',1234);
A = r.randn(6,6);

n = length(A);
p = [1:n]';
M = A;

for i = 1:n-1
   for m = i:n
      for k = m+1:n
          if abs(M(k,i)) >= abs(M(m,i))
             blah = M(i,:); M(i,:) = M(k,:); M(k,:) = blah;
             blah = p(i); p(i) = p(k); p(k) = blah;
          end
              
       end
   end
   
   for j = i+1:n
      a = M(j,i)/M(i,i);
      M(j,i:end) = M(j,i:end) - a*M(i,i:end);
      M(j,i) = a;
   end
end

A1 = p; A2 = M;

b = [1;0;1;0;1;0];
L = tril(M) - diag(diag(M)) + eye(length(M));
U = triu(M);

z = b(p); A3 = z;
y = Forsub([L,z]); A4 = y;
x = Backsub([U,y]); A5 = x;

A = [1.1,0.2,-0.2,0.5;
     0.2,0.9,0.5,0.3;
     0.1,0.0,1.0,0.4;
     0.1,0.1,0.1,1.2];
b = [1;0;1;0];

y = zeros(length(b),1);
M = eye(length(b)) - A;

for i = 1:100000
   y = M*y + b;
   if max(abs(A*y - b)) < 1e-2
       A6 = i; A7 = max(abs(A*y - b));
       break
   end
end

y = zeros(length(b),1);
M = eye(length(b)) - A;

for i = 1:100000
   y = M*y + b;
   if max(abs(A*y - b)) < 1e-4
       A8 = i; A9 = max(abs(A*y - b));
       break
   end
end

y = zeros(length(b),1);
M = eye(length(b)) - A;

for i = 1:100000
   y = M*y + b;
   if max(abs(A*y - b)) < 1e-6
       A10 = i; A11 = max(abs(A*y - b));
       break
   end
end

save('CP1Sols.mat','A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11')

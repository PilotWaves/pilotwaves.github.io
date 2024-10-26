%%%  Problem 1 %%%
r = RandStream('mt19937ar','Seed',1234);
A = r.randn(10,6);

Q = eye(length(A));  R = A;

for i = 1:length(A(1,:))
   w = R(i:end,i);
   w(1) = w(1) + sign(w(1))*norm(w);
   w = w/norm(w);
   R(i:end,i:end) = R(i:end,i:end) - 2*w*(w'*R(i:end,i:end));
   Q(:,i:end) = Q(:,i:end) - 2*(Q(:,i:end)*w)*w';
end

A1 = Q; A2 = R;


%%%  Problem 2 %%%

r = RandStream('mt19937ar','Seed',1234);
A = r.randn(10,6);
b = r.randn(10,1);

A3 = A'*A; A4 = A'*b; A5 = A3\A4; A6 = abs(b-A*A5);


%%%  Problem 3 %%%

r = RandStream('mt19937ar','Seed',1234);
A = r.randn(10,10)+1i*r.randn(10,10);

[V,D] = eig(A);

DiagD = diag(D);
RealDiag = real(DiagD);
p = [1:length(RealDiag)];

for i = 1:length(RealDiag)
   for j = i+1:length(RealDiag)
       if RealDiag(j) < RealDiag(i)
          blah = p(i); p(i) = p(j); p(j) = blah;
          blah = RealDiag(i); RealDiag(i) = RealDiag(j); RealDiag(j) = blah;
       end
   end
end

DiagD = DiagD(p); D = diag(DiagD); V = V(:,p);
A7 = real(V); A8 = imag(V); A9 = real(D); A10 = imag(D);

save('CP3Sols.mat','A1','A2','A3','A4','A5','A6','A7','A8','A9','A10');
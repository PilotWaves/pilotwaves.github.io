function [y] = LinearJacobi(y,A,b)

L = tril(A,-1);
U = triu(A,+1);
D = diag(diag(A));

c = inv(D)*b;

y = -inv(D)*(L + U)*y + c;
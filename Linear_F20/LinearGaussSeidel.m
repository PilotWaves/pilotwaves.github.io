function [y] = LinearGaussSeidel(y,A,b)

L = tril(A,-1);
U = triu(A,+1);
D = diag(diag(A));

y = -inv(L+D)* U*y + inv(L+D)*b;
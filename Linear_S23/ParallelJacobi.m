function [y] = ParallelJacobi(y,A,b)

L = tril(A,-1);
U = triu(A,+1);
D = diag(A);  % this is just a vector

FirstTermMultiplication = -(L + U)*y;

parfor i = 1:length(b)
    y(i) = FirstTermMultiplication(i)/D(i) + b(i)/D(i);
end
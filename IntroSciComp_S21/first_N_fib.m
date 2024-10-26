function fib = first_N_fib(F1, F2, N)

fib = zeros(1,N);
fib(1) = F1;
fib(2) = F2;

for n = 3:N
    fib(n) = fib(n-1) + fib(n-2);
end

end
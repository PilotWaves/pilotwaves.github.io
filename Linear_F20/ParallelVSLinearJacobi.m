ParallelTime = 0*[0:4];
LinearTime = ParallelTime;

for n = 0:length(ParallelTime)-1

A = diag([1:10^n]) + .01*randn(10^n);
b = ones(length(A),1);

y= zeros(length(b),1);

E = 1e-8;
T = 0;

tic
while max(abs(A*y - b)) > E
    y = ParallelJacobi(y,A,b);
end
ParallelTime(n+1) = toc;


y= zeros(length(b),1);

tic
while max(abs(A*y - b)) > E
    y = LinearJacobi(y,A,b);
end
LinearTime(n+1) = toc;

end

set(0,'defaultAxesFontSize',20)
semilogy([0:length(ParallelTime)-1],ParallelTime,[0:length(ParallelTime)-1],LinearTime,'Linewidth',5)
xlabel('n (powers of 10)')
ylabel('Runtime (log scale)')

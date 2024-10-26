gamma = 0.44;

N = round(1/gamma) + 1;

v = zeros(1,N);

del = 31/2;
CK = 1/7;

v(1) = 2*pi*5/(3*del);
s = 2*pi/del;

set(0,'DefaultAxesFontSize',14)

for m = 1:10
for n = 1:N-1
    
    v(n+1) = v(n);

        if gamma*n < 1
            v(n+1) = v(n) + CK*exp(1 - 1/(1-(n*gamma)^2))*exp(-m*s*del/(8.4*pi));
        end
end

% set(gca,'Yticklabel',[]) 
% set(gca,'Xticklabel',[])
plot(m*3,v(n+1)/v(1),'p','color',[0 0.7 0],'MarkerSize',20)
axis([0 30 1 1.15])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
hold on

end

hold off
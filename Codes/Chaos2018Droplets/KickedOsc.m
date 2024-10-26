gamma = 0.44;

N = round(1/gamma) + 1;

v = zeros(1,N);

del = 31/2;
CK = 1/7;

v(1) = 2*pi*5/(3*del);
s = (8/6)*2*pi/del;

set(0,'DefaultAxesFontSize',14)

for M = 1:11
for n = 1:N-1
    
    v(n+1) = v(n);

    if gamma*n < 1
    
        for m = 1:M
            v(n+1) = v(n+1) + CK*exp(1 - 1/(1-(n*gamma)^2))*exp(-m*s*del/(4.2*pi));
        end
        
    end
    
end

v(N)/v(1)

% set(gca,'Yticklabel',[]) 
% set(gca,'Xticklabel',[])
plot(m+1,v(N)/v(1),'p','color',[0 0.7 0],'MarkerSize',20)
axis([0 12 1 1.25])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
hold on

end

% set(gca,'Yticklabel',[]) 
% set(gca,'Xticklabel',[])
plot(1,1,'p','color',[0 0.7 0],'MarkerSize',20)
axis([0 12 1 1.25])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
hold off

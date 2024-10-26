N = 100;
t = [1:N];

v = zeros(1,N);
theta = v;

v(1) = 1;
s = 1;

gamma = 0.71;

set(0,'DefaultAxesFontSize',14)

for M = 1:11
for n = 1:N-1
    
    v(n+1) = v(n);

    if gamma*n < 1
    
        for m = 1:M
            v(n+1) = v(n+1) + (1/gamma)*exp(-1/(1-(n*gamma)^2))*exp(-gamma*m*s)*cos((1-gamma)*theta(n));
        end
        
    end
    
    %v(n+1) = mod(v(n+1),2*pi);
    theta(n+1) = mod(theta(n) + v(n+1),2*pi);
    
end

% set(gca,'Yticklabel',[]) 
% set(gca,'Xticklabel',[])
plot(m+1,v(N),'p','color',[0 0.7 0],'MarkerSize',20)
%axis([0 12 1 1.25])
hold on

end

% set(gca,'Yticklabel',[]) 
% set(gca,'Xticklabel',[])
plot(1,1,'p','color',[0 0.7 0],'MarkerSize',20)
%axis([0 12 1 1.25])
hold off

% figure(1)
% h = plot(t,v,'-');
% set(h(1),'linewidth',1)
% xlabel('n')
% ylabel('v')

% figure(2)
% g = plot(t,mod(theta,2*pi),'--');
% set(g(1),'linewidth',1)
% xlabel('n')
% ylabel('\theta^i - \theta^j')

clear all

mu = 3.2;
nu = 2*mu/(1+mu^2);

x1 = [-nu*mu:.001:-nu];
x2 = [-nu:.001:nu];
x3 = [nu:.001:nu*mu];

y1 = (1+mu)*x1/(1-mu) + (nu + mu*nu)/(1-mu) - mu*nu;
y2 = mu*x2;
y3 = (1+mu)*x3/(1-mu) - (nu + mu*nu)/(1-mu) + mu*nu;

x = [x1 x2 x3];

set(0,'defaultlinelinewidth',3)
set(0,'DefaultAxesFontSize',24)

h = plot(x,x,'r--',x1,y1,'b',x2,y2,'b',x3,y3,'b');
set(h(1),'linewidth',2) %y = x
xlabel('x');
ylabel('y');
xL = xlim;
yL = ylim;
line(xL, [0 0],'color','k','linewidth',1) %x-axis 
line([0 0], yL,'color','k','linewidth',1) %y-axis
g = legend([h(1) h(2)],{'y=x','y_f'},'Location','northwest');
legend boxoff
set(g,'FontSize',20,'FontWeight','bold')
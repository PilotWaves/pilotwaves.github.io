clear all

mu = 2;
nu = 1/(1+mu);

x1 = [nu-1:.001:nu-1/2];
x2 = [nu-1/2:.001:nu];

y1 = -1 + nu + mu*(1+x1-nu);
y2 = -1 + nu + mu*(nu-x2);

x = [x1 x2];

set(0,'defaultlinelinewidth',3)
set(0,'DefaultAxesFontSize',24)

h = plot(x,x,'r--',x1,y1,'b',x2,y2,'b');
set(h(1),'linewidth',2) %y = x
xlabel('x');
ylabel('y')
xL = xlim;
yL = ylim;
line(xL, [0 0],'color','k','linewidth',1) %x-axis 
line([0 0], yL,'color','k','linewidth',1) %y-axis
g = legend([h(1) h(2)],{'y=x','y_g-1'},'Location','northwest');
legend boxoff
set(g,'FontSize',20,'FontWeight','bold')
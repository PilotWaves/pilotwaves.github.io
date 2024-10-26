theta_1 = pi/2; %initial condition of the first
theta_2 = pi/2; %initial condition of the second
duration = 10; %duration in seconds

L1 = 2; %length of the first in meters
L2 = 1; %length of the second in meters

[T,Y] = ode45(@DoublePendulum,[0 duration],[theta_1 theta_2 0 0]);

y1 = -L1*cos(Y(:,1));
y2 = -(L1*cos(Y(:,1))+L2*cos(Y(:,2)));
x1 = L1*sin(Y(:,1));
x2 = L1*sin(Y(:,1))+L2*sin(Y(:,2));

for t = 1:length(y1)
h = plot([0 x1(t)],[0 y1(t)],'k',[x1(t) x2(t)],[y1(t) y2(t)],'k',x1(t),y1(t),'.',0,0,'.',x2(t),y2(t),'.',x2(1:t),y2(1:t),'r');
set(h(1),'linewidth',2);
set(h(2),'linewidth',2);
set(h(3),'MarkerSize',20);
set(h(4),'MarkerSize',30);
set(h(5),'MarkerSize',20);
title(['t=' num2str(duration*t/(length(y1)))]);
axis([-3 3 -3 3])
pause;
end

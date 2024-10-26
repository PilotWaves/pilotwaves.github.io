N = 10^7;
t = [1:N];

v = zeros(1,N);
theta = v;

theta1 = [0:2*pi/100:2*pi];
x1 = 1.25*cos(theta1);
y1 = 1.25*sin(theta1);
x2 = 0.75*cos(theta1);
y2 = 0.75*sin(theta1);

set(0,'DefaultAxesFontSize',14)

%gamma = .00001;
del = 31/2;
K = -(pi*exp(del^2/8.4))/(sin(pi*del));
mu = 1/K;

v(1) = 5/(3*del);

% writerObj = VideoWriter('ConstantSpeed.avi');
% writerObj.FrameRate = 10;
% open(writerObj);

for n = 1:N-1

    v(n+1) = mu*(v(n) + K*sin(del*v(n))*exp(-v(n)^2*del^2/(8.4*pi^2)));
    theta(n+1) = theta(n) + v(n+1);
    
%     x = cos(theta(1:n));
%     y = sin(theta(1:n));
    
%     h = plot(x,y,'o',x1,y1,'k',x2,y2,'k',x(n),y(n),'r.');
%      h = plot(x(n),y(n),'.',x1,y1,'k',x2,y2,'k');
%     set(h(1),'MarkerSize',80);
%     set(h(2),'LineWidth',4);
%     set(h(3),'LineWidth',4);
%    set(h(4),'MarkerSize',20);
%     axis([-1.5 1.5 -1.5 1.5],'square');
%     axis off
     
%    txt = ['   ',num2str(n)];
%     text(x(n),y(n),txt,'FontSize',18)
%     hold on
    

%     frame = getframe(gcf);
%     writeVideo(writerObj,frame);
    
end

%     x = cos(theta);
%     y = sin(theta);

%    txt = ['   ',num2str(6)];
%     text(x(6),y(6),txt,'FontSize',18)
% 
% hold off

% close(writerObj)
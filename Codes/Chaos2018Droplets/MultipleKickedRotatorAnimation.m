N = 100;
M = 5;
t = [1:N];

v = zeros(M,N);
theta = v;

theta1 = [0:2*pi/100:2*pi];
x1 = 1.25*cos(theta1);
y1 = 1.25*sin(theta1);
x2 = 0.75*cos(theta1);
y2 = 0.75*sin(theta1);

set(0,'DefaultAxesFontSize',14)

del = 31/2;
K = -(pi*exp(del^2/8.4))/(sin(pi*del));
mu = 1/(3*K);

v(:,1) = 5/(3*del);
theta(:,1) = [0:2*pi/M:2*pi - 2*pi/M];

writerObj = VideoWriter('MultipleSpeed.avi');
writerObj.FrameRate = 5;
open(writerObj);

set(0,'DefaultAxesColorOrder',[[0 0 1]; [0 1 0]; [1 0 0]; [0 0 0]; [1 1 0]])

for n = 1:N-1
for i = 1:M
        
    sum = 0;
    for m = 1:M
       if m ~= i
           sum = sign(mod(theta(i,n),2*pi) - mod(theta(m,n),2*pi))*exp(-(del^2/(8.4*pi^2))*(mod(theta(i,n) - theta(m,n),2*pi))^2);
       end
    end

    v(i,n+1) = mu*(v(i,n) + K*sin(del*v(i,n))*exp(-v(i,n)^2*del^2/(8.4*pi^2)) + K*sum);
    theta(i,n+1) = theta(i,n) + v(i,n+1);
    
    x = cos(theta(i,n+1));
    y = sin(theta(i,n+1));

     h = plot(x,y,'.',x1,y1,'k',x2,y2,'k');
    set(h(1),'MarkerSize',20);
    set(h(2),'LineWidth',4);
    set(h(3),'LineWidth',4);
    axis([-1.5 1.5 -1.5 1.5],'square');
    axis off
    hold on
     
%    txt = ['    ',num2str(n)];
%     text(x(n),y(n),txt)
%     hold on
    
end

    frame = getframe(gcf);
    writeVideo(writerObj,frame);
    hold off

end

%     x = cos(theta);
%     y = sin(theta);

%    txt = ['    ',num2str(6)];
%     text(x(6),y(6),txt)
% 
hold off

 close(writerObj)
clear all

N = 10000;

y1 = zeros(1,N+1);
x1 = zeros(1,N+1);

x2 = x1;
y2 = y1;

 writerObj = VideoWriter('GlobalBifurcation_Origin.avi');
 writerObj.FrameRate = 5;
 open(writerObj);

set(0,'DefaultAxesFontSize',14)

for sigma = 0.4:.01:0.95

y1(1) = 0;
x1(1) = 0.5;

y2(1) = 0;
x2(1) = 1.5;

for n = 1:N

x1(n+1) = x1(n) - 4*sigma*(1-x1(n))*y1(n);
y1(n+1) = 0.2*(y1(n)+x1(n)*(2-x1(n)));

x2(n+1) = x2(n) - 4*sigma*(1-x2(n))*y2(n);
y2(n+1) = 0.2*(y2(n)+x2(n)*(2-x2(n)));

end

h = plot(x1,y1,'r.',x2,y2,'.');
set(h(1),'MarkerSize',1);
set(h(2),'MarkerSize',1);
axis([-1.2 3.2 -0.75 0.26])
alw = 0.75;    % AxesLineWidth
fsz = 14;      % Fontsize
xlabel('x')
ylabel('y')
title(['\sigma = ' num2str(sigma)]);
%pause

 frame = getframe(gcf);
 writeVideo(writerObj,frame);
end

close(writerObj);
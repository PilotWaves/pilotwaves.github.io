clear all

N = 10000;

y1 = zeros(1,N+1);
x1 = zeros(1,N+1);

x2 = x1;
y2 = y1;

%  writerObj = VideoWriter('GlobalBifurcation1_Origin.avi');
%  writerObj.FrameRate = 5;
%  open(writerObj);

set(0,'DefaultAxesFontSize',14)

for sigma = 0:.0001:0.025

y1(1) = 0.1;
x1(1) = 0.5;

y2(1) = 0.1;
x2(1) = 1.5;

for n = 1:N

x1(n+1) = x1(n) - sigma*(216*(x1(n)^2 - 9*x1(n)/4 + 5/4)/11)*y1(n);
y1(n+1) = 0.1*(y1(n)+216*(x1(n)^3/3 - 9*x1(n)^2/8 + 5*x1(n)/4)/11);

x2(n+1) = x2(n) - sigma*(216*(x2(n)^2 - 9*x2(n)/4 + 5/4)/11)*y2(n);
y2(n+1) = 0.1*(y2(n)+216*(x2(n)^3/3 - 9*x2(n)^2/8 + 5*x2(n)/4)/11);

end

h = plot(x1,y1,'r.',x2,y2,'.');
set(h(1),'MarkerSize',4);
set(h(2),'MarkerSize',4);
%axis([-1.2 3.2 -0.75 0.26])
alw = 0.75;    % AxesLineWidth
fsz = 14;      % Fontsize
xlabel('x')
ylabel('y')
title(['\lambda = ' num2str(sigma)]);
pause

%  frame = getframe(gcf);
%  writeVideo(writerObj,frame);
end

%close(writerObj);
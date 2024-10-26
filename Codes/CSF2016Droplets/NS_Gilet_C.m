clear all

N = 100000;

w1 = zeros(1,N+1);
x1 = zeros(1,N+1);

x2 = x1;
w2 = w1;

B = pi/3;

mu = 0.5;

 %writerObj = VideoWriter('GlobalBifurcation_Origin.avi');
 %writerObj.FrameRate = 5;
 %open(writerObj);



for C = 0.45:.01:0.9

w1(1) = 0;
x1(1) = -1.6;

w2(1) = 0;
x2(1) = -1.55;

for n = 1:N

psi1 = (cos(B)*sin(3*x1(n))+sin(B)*sin(5*x1(n)))/sqrt(pi);
dpsi1 = (3*cos(B)*cos(3*x1(n))+5*sin(B)*cos(5*x1(n)))/sqrt(pi);

w1(n+1) = mu*(w1(n)+psi1);
x1(n+1) = x1(n) - C*w1(n)*dpsi1;

psi2 = (cos(B)*sin(3*x2(n))+sin(B)*sin(5*x2(n)))/sqrt(pi);
dpsi2 = (3*cos(B)*cos(3*x2(n))+5*sin(B)*cos(5*x2(n)))/sqrt(pi);

w2(n+1) = mu*(w2(n)+psi2);
x2(n+1) = x2(n) - C*w2(n)*dpsi2;

end

h = plot(x1,w1,'r.',x2,w2,'.');
set(h(1),'MarkerSize',1);
set(h(2),'MarkerSize',1);
axis([-2.3 -0.8 -0.3 0.4])
%axis([-0.35 0.35 -0.3 0.4])
alw = 0.75;    % AxesLineWidth
fsz = 14;      % Fontsize
xlabel('x')
ylabel('w')
title(['C = ' num2str(C)]);
pause

 %frame = getframe(gcf);
 %writeVideo(writerObj,frame);
end

%close(writerObj);
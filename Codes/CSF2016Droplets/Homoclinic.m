clear all

N = 10000;

w1 = zeros(1,N+1);
x1 = zeros(1,N+1);
x2 = x1;
w2 = w1;

B = pi/3;
%mu = .908312496093;
%mu = .915;
C = 0.05;

for mu = 0.913:.001:.915
xc = pi/2;

psic = (cos(B)*sin(3*xc)+sin(B)*sin(5*xc))/sqrt(pi);

wc = mu*psic/(1-mu);

w1(1) = wc-.5;
x1(1) = xc-.1;
w2(1) = wc-.5;
x2(1) = xc+.1;

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

h = plot(x1,w1,'.',[xc-.2 xc+.2],[wc wc],'m--',[xc xc],[wc-1.5 wc+1.5],'g--',x2,w2,'r.');
set(h(1),'MarkerSize',4);
set(h(2),'linewidth',1);
set(h(3),'linewidth',1);
set(h(1),'MarkerSize',4);
%axis([xb-.1 xc+.1 wb-1 wc+1])
alw = 0.75;    % AxesLineWidth
fsz = 12;      % Fontsize
xlabel('x')
ylabel('w')
title(['mu = ' num2str(mu)]);
pause;
end

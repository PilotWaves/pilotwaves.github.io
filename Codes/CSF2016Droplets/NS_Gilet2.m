clear all

N = 500;
M = 203;

w1 = zeros(1,N+1);
x1 = zeros(1,N+1);

w2 = w1;
x2 = x1;

w3 = zeros(1,M+1);
x3 = zeros(1,M+1);

B = 5*pi/6;
mu = .999;
C = 0.05;

w1(1) = 0;
x1(1) = 0.2;

w2(1) = 0;
x2(1) = 0.3;

w3(1) = 0;
x3(1) = 0.2644;

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

for n = 1:M
psi3 = (cos(B)*sin(3*x3(n))+sin(B)*sin(5*x3(n)))/sqrt(pi);
dpsi3 = (3*cos(B)*cos(3*x3(n))+5*sin(B)*cos(5*x3(n)))/sqrt(pi);

w3(n+1) = mu*(w3(n)+psi3);
x3(n+1) = x3(n) - C*w3(n)*dpsi3; 
end

h = plot(x1,w1,'.',x2,w2,'r.',x1(1),w1(1),'b*',x2(1),w2(1),'r*',x3,w3,'k--');
set(h(1),'MarkerSize',8);
set(h(2),'MarkerSize',8);
set(h(3),'MarkerSize',10);
set(h(4),'MarkerSize',10);
set(h(5),'linewidth',1)
alw = 0.75;    % AxesLineWidth
fsz = 12;      % Fontsize
xlabel('x')
ylabel('w')
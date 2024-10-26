clear all

N = 10000;

w1 = zeros(1,N+1);
x1 = w1;

w2 = w1;
x2 = x1;

w3 = w1;
x3 = x1;

w4 = w1;
x4 = x1;

B = pi/3;
mu = .99;
C = 0.05;

w1(1) = 0;
x1(1) = pi + .1;

w2(1) = 0.1;
x2(1) = pi;

w3(1) = 0;
x3(1) = pi - .1;

w4(1) = -.1;
x4(1) = pi;

for n = 1:N

psi1 = (cos(B)*sin(3*(x1(n)-pi/2))+sin(B)*sin(5*(x1(n)-pi/2)))/sqrt(pi);
dpsi1 = (3*cos(B)*cos(3*(x1(n)-pi/2))+5*sin(B)*cos(5*(x1(n)-pi/2)))/sqrt(pi);

w1(n+1) = mu*(w1(n)+psi1);
x1(n+1) = mod(x1(n) - C*w1(n)*dpsi1,2*pi);

psi2 = (cos(B)*sin(3*(x2(n)-pi/2))+sin(B)*sin(5*(x2(n)-pi/2)))/sqrt(pi);
dpsi2 = (3*cos(B)*cos(3*(x2(n)-pi/2))+5*sin(B)*cos(5*(x2(n)-pi/2)))/sqrt(pi);

w2(n+1) = mu*(w2(n)+psi2);
x2(n+1) = mod(x2(n) - C*w2(n)*dpsi2,2*pi);

psi3 = (cos(B)*sin(3*(x3(n)-pi/2))+sin(B)*sin(5*(x3(n)-pi/2)))/sqrt(pi);
dpsi3 = (3*cos(B)*cos(3*(x3(n)-pi/2))+5*sin(B)*cos(5*(x3(n)-pi/2)))/sqrt(pi);

w3(n+1) = mu*(w3(n)+psi3);
x3(n+1) = mod(x3(n) - C*w3(n)*dpsi3,2*pi);

psi4 = (cos(B)*sin(3*(x4(n)-pi/2))+sin(B)*sin(5*(x4(n)-pi/2)))/sqrt(pi);
dpsi4 = (3*cos(B)*cos(3*(x4(n)-pi/2))+5*sin(B)*cos(5*(x4(n)-pi/2)))/sqrt(pi);

w4(n+1) = mu*(w4(n)+psi4);
x4(n+1) = mod(x4(n) - C*w4(n)*dpsi4,2*pi);

end

x = [x1 x2 x3 x4];
hist(x,4*N)
alw = 0.75;    % AxesLineWidth
fsz = 28;      % Fontsize
xlabel('x')
ylabel('Number of Occurances')

% h = plot(x1,w1,'b.',x2,w2,'b.',x3,w3,'b.',x4,w4,'b.');
% set(h(1),'MarkerSize',1);
% set(h(2),'MarkerSize',1);
% set(h(3),'MarkerSize',1);
% set(h(4),'MarkerSize',1);
% alw = 0.75;    % AxesLineWidth
% fsz = 12;      % Fontsize
% xlabel('x')
% ylabel('w')
% axis([0 2*pi -6 6])
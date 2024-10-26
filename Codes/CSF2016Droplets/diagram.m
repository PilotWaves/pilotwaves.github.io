clear all

N = 100;

w = zeros(1,N+1);
x = w;

B = pi/3;
mu = .7;
C = 0.05;

w(1) = 0;
x(1) = .35;

X = [x(1)-pi/2:pi/400:x(1)+pi/2];
W = 0*X;

psi = (cos(B)*sin(3*x(1))+sin(B)*sin(5*x(1)))/sqrt(pi);

for n = 1:N-1

f = plot(X,W,'b',x(n),w(n)*psi+5/9,'bo',x(n),w(n)*psi+5/9,'b.',[x(n) x(n)],...
    [0 w(n)*psi],'k');
set(f(1),'linewidth',4,'LineSmoothing','on');
set(f(2),'MarkerSize',80);
set(f(3),'MarkerSize',20);
set(f(4),'linewidth',2);
warning('off','last')
%light('Position',[x(n+1),w(n+1)*psi+5/6,0])
axis([x(n)-pi/8 x(n)+pi/8 -1 3])
hold on

dpsi= (3*cos(B)*cos(3*x(n))+5*sin(B)*cos(5*x(n)))/sqrt(pi);

w(n+1) = mu*(w(n)+psi);
x(n+1) = x(n) - C*w(n)*dpsi;

AW = (w(n) + psi)*(cos(B)*sin(3*(X))+sin(B)*sin(5*(X)))/sqrt(pi);

g = plot(X,AW,'g--',x(n),(w(n)+psi)*psi+3/5,'go',x(n),(w(n)+psi)*psi+3/5,'g.',...
    [x(n) x(n)],[0 (w(n)+psi)*psi],'k');
set(g(1),'linewidth',4);
set(g(2),'MarkerSize',80);
set(g(3),'MarkerSize',20);
set(g(4),'linewidth',2);
axis([x(n)-pi/8 x(n)+pi/8 -1 3])
hold on

X = [x(n+1)-pi/2:pi/400:x(n+1)+pi/2];
W = w(n+1)*(cos(B)*sin(3*(X))+sin(B)*sin(5*(X)))/sqrt(pi);
psi = (cos(B)*sin(3*x(n+1))+sin(B)*sin(5*x(n+1)))/sqrt(pi);

h = plot(X,W,'r',x(n+1),w(n+1)*psi+3/5,'ro',x(n+1),w(n+1)*psi+3/5,...
    'r.',[x(n+1) x(n+1)],[0 w(n+1)*psi],'k',[x(n)-pi/8 x(n)+pi/8],[0 0],'k--');
set(h(1),'linewidth',4,'LineSmoothing','on');
set(h(2),'MarkerSize',80);
set(h(3),'MarkerSize',20);
set(h(4),'linewidth',2);
legend([f(3),g(3),h(3),f(1),g(1),h(1)],'Just before impact n','Just after impact n',...
    'Just before impact n+1','w_n\psi(x)','[w_n+\psi(x_n)]\psi(x)','w_{n+1}\psi(x)')
warning('off','last')
%light('Position',[x(n+1),w(n+1)*psi+5/6,0])
axis([x(n)-pi/8 x(n)+pi/8 -1 3])
hold off
pause;

end
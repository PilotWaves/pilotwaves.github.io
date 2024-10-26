clear all

N = 10000;

w = zeros(1,N+1);
x = w;

B = pi/3;
mu = .94;
C = 0.05;

w(1) = 0;
x(1) = 1.57079632;
%x(1) = .35;

X = [x(1)-pi/4:pi/400:x(1)+pi/4];
W = 0*X;

M = 100;
strat = [-20:20/M:0]';
blah = strat*ones(1,length(X))+ones(length(strat),1)*W;

h = plot(X,W,x(1),w(1)+5/6,'bo');
set(h(1),'linewidth',2,'LineSmoothing','on');
set(h(2),'MarkerSize',20,'MarkerFaceColor',[28/255,169/255,201/255]);
hold on

for m = 1:M
fill([X,fliplr(X)],[blah(m,:),fliplr(blah(m+1,:))],[28/255,169/255,201/255],...
    'edgecolor','none','facealpha',1-(.4*m/M),'LineSmoothing','on');
hold on
end
warning('off','last')
axis([x(1)-pi/4 x(1)+pi/4 -10 10])
hold off

psi = (cos(B)*sin(3*x(1))+sin(B)*sin(5*x(1)))/sqrt(pi);

for n = 1:N-1

dpsi= (3*cos(B)*cos(3*x(n))+5*sin(B)*cos(5*x(n)))/sqrt(pi);

w(n+1) = mu*(w(n)+psi);
x(n+1) = x(n) - C*w(n)*dpsi;

AW = (w(n) + psi)*(cos(B)*sin(3*(X))+sin(B)*sin(5*(X)))/sqrt(pi);
blah = strat*ones(1,length(X))+ones(length(strat),1)*AW;

h = plot(X,AW,x(n),(w(n)+psi)*psi+5/6,'bo');
set(h(1),'linewidth',2,'LineSmoothing','on');
set(h(2),'MarkerSize',20,'MarkerFaceColor',[28/255,169/255,201/255]);
hold on
for m = 1:M
fill([X,fliplr(X)],[blah(m,:),fliplr(blah(m+1,:))],[28/255,169/255,201/255],...
    'edgecolor','none','facealpha',1-(.4*m/M),'LineSmoothing','on');
hold on
end
warning('off','last')
%light('Position',[x(n),(w(n)+psi)*psi+5/6,0])
axis([x(n)-pi/4 x(n)+pi/4 -10 10])
hold off

Ix = (x(n+1)+x(n))/2;
Ix0 = (Ix + x(n))/2;
Ix1 = (Ix + x(n+1))/2;
X = [x(n+1)-pi/4:pi/400:x(n+1)+pi/4];
IX = [Ix-pi/4:pi/400:Ix+pi/4];
IX0 = [Ix0-pi/4:pi/400:Ix0+pi/4];
IX1 = [Ix1-pi/4:pi/400:Ix1+pi/4];
W = w(n+1)*(cos(B)*sin(3*(X))+sin(B)*sin(5*(X)))/sqrt(pi);
IW = (w(n+1)+w(n)+psi)*(cos(B)*sin(3*(IX))+sin(B)*sin(5*(IX)))/(2*sqrt(pi));
IW0 = (IW(101)+w(n)+psi)*(cos(B)*sin(3*(IX0))+sin(B)*sin(5*(IX0)))/(2*sqrt(pi));
IW1 = (w(n+1)+IW(101))*(cos(B)*sin(3*(IX1))+sin(B)*sin(5*(IX1)))/(2*sqrt(pi));
psi = (cos(B)*sin(3*x(n+1))+sin(B)*sin(5*x(n+1)))/sqrt(pi);

blah = strat*ones(1,length(X))+ones(length(strat),1)*IW0;

h = plot(IX0,IW0,Ix0,mu*(w(n)*psi+3),'bo');
set(h(1),'LineSmoothing','on');
set(h(2),'MarkerSize',20,'MarkerFaceColor',[20/255,160/255,201/255]);
hold on
for m = 1:M
fill([IX0,fliplr(IX0)],[blah(m,:),fliplr(blah(m+1,:))],[28/255,169/255,201/255],...
    'edgecolor','none','facealpha',1-(.45*m/M),'LineSmoothing','on');
hold on
end
warning('off','last')
%light('Position',[Ix0,mu*(w(n)*psi+4),0])
axis([Ix0-pi/4 Ix0+pi/4 -10 10])
hold off

blah = strat*ones(1,length(X))+ones(length(strat),1)*IW;

h = plot(IX,IW,Ix,mu*(w(n+1)*psi+4),'bo');
set(h(1),'LineSmoothing','on');
set(h(2),'MarkerSize',20,'MarkerFaceColor',[20/255,160/255,201/255]);
hold on
for m = 1:M
fill([IX,fliplr(IX)],[blah(m,:),fliplr(blah(m+1,:))],[28/255,169/255,201/255],...
    'edgecolor','none','facealpha',1-(.45*m/M),'LineSmoothing','on');
hold on
end
warning('off','last')
%light('Position',[Ix,mu*(w(n+1)*psi+4),0])
axis([Ix-pi/4 Ix+pi/4 -10 10])
hold off

if mod(n,100) == 0
pause
end

blah = strat*ones(1,length(X))+ones(length(strat),1)*IW1;

h = plot(IX1,IW1,Ix1,mu*(w(n+1)*psi+3),'bo');
set(h(1),'LineSmoothing','on');
set(h(2),'MarkerSize',20,'MarkerFaceColor',[20/255,160/255,201/255]);
hold on
for m = 1:M
fill([IX1,fliplr(IX1)],[blah(m,:),fliplr(blah(m+1,:))],[28/255,169/255,201/255],...
    'edgecolor','none','facealpha',1-(.45*m/M),'LineSmoothing','on');
hold on
end
warning('off','last')
%light('Position',[Ix1,mu*(w(n+1)*psi+2),0])
axis([Ix1-pi/4 Ix1+pi/4 -10 10])
hold off

blah = strat*ones(1,length(X))+ones(length(strat),1)*W;

h = plot(X,W,x(n+1),w(n+1)*psi+5/6,'bo');
set(h(1),'linewidth',2,'LineSmoothing','on');
set(h(2),'MarkerSize',20,'MarkerFaceColor',[28/255,169/255,201/255]);
hold on
for m = 1:M
fill([X,fliplr(X)],[blah(m,:),fliplr(blah(m+1,:))],[28/255,169/255,201/255],...
    'edgecolor','none','facealpha',1-(.4*m/M),'LineSmoothing','on');
hold on
end
warning('off','last')
%light('Position',[x(n+1),w(n+1)*psi+5/6,0])
axis([x(n+1)-pi/4 x(n+1)+pi/4 -10 10])
hold off

end
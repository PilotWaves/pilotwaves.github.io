clear

%%%External Parameters%%%
T1 = 20;  %period of I1 and I4
T2 = 40;  %period of I2 and I3
Total = 80;  %Total time of clock
freq = 1/9;  %internal frequency
%%%%%%%%%%%%%%%%%%%%%%%%%
N = Total/freq;
t = 0:freq:Total;

%%%Internal Parameters%%%
mu = 3.2;
nu = 2*mu/(1+mu^2);
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Initialization%%%
f = zeros(1,N+1);  %VT1
g = f;  %VT2
h1 = f;  %U1
h2 = f;  %U2
%%%%%%%%%%%%%%%%%%%%

%%%Inputs%%%
I1 = f;
I1(T1/(2*freq)+1:T1/freq+1) = 1;% - rand*1e-8;
I1(3*T1/(2*freq)+1:2*T1/freq+1) = 1;% - rand*1e-8;
I1(5*T1/(2*freq)+1:3*T1/freq+1) = 1;% - rand*1e-8;
I1(7*T1/(2*freq)+1:4*T1/freq+1) = 1;% - rand*1e-8;
I2 = f;
I2(T2/(2*freq)+2:T2/freq+1) = 1 + 1e-9;% - rand*1e-8;
I2(3*T2/(2*freq)+2:2*T2/freq+1) = 1 - 7e-10;% - rand*1e-8;
I3 = I2;
I4 = I1;
h1(1) = 1;
h2(1) = 1;
f(1) = 0;
g(1) = -1;
%%%%%%%%%%%%

%%%Test%%%
% subplot(2,1,1)
% plot(t,I1)
% subplot(2,1,2)
% plot(t,I2)
%%%%%%%%%%

for n = 1:N
    
    %%%Adding noise into inputs%%%
%         I1(n+1) = I1(n+1) + rand*1e-8;
%         I2(n+1) = I2(n+1) + rand*1e-8;
%         I3(n+1) = I3(n+1) + rand*1e-8;
%         I4(n+1) = I4(n+1) + rand*1e-8;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   
    %%%VT1 calculations%%%
    
    if  f(n) >= -nu*mu && f(n) < -nu
        
        yf = (1+mu)*f(n)/(1-mu) + nu*(1+mu)/(1-mu) - mu*nu;
        
    else if f(n) >= -nu && f(n) <= nu
            
            yf = mu*f(n);
            
        else
            
            yf = (1+mu)*f(n)/(1-mu) - nu*(1+mu)/(1-mu) + mu*nu;
            
        end
        
    end
    
    f(n+1) = abs(I1(n+1) - I2(n+1)) + I1(n+1)*I2(n+1)*yf;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%VT2 calculations%%%
    
    if abs(g(n)-1/3) <= 1/2
    
        yg = 1/3 + 2*abs(g(n)-1/3);
        
    else
        
        yg = 1/3 + 2*(1-abs(g(n)-1/3));
        
    end
    
    g(n+1) = abs(I3(n+1) - I4(n+1)) - 1 + I3(n+1)*I4(n+1)*yg;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%NOR Output calculations%%%
    
h1(n+1) = (1-I1(n+1))*(1-I2(n+1));
h2(n+1) = (1-I3(n+1))*(1-I4(n+1));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

set(0,'defaultlinelinewidth',2)
set(0,'DefaultAxesFontSize',24)

subplot(4,1,1);
plot(t,I1);
axis([0 80 0 1.08]);
ylabel('I1')
subplot(4,1,2);
plot(t,I2);
axis([0 80 0 1.08]);
ylabel('I2')
subplot(4,1,3);
plot(t,f);
axis([0 80 -2.16 2.16]);
ylabel('VT1')
subplot(4,1,4);
plot(t,h1);
axis([0 80 -1.08 1.08]);
ylabel('U1')
xlabel('Time in ms')

figure

subplot(4,1,1);
plot(t,I3);
axis([0 80 0 1.08]);
ylabel('I3')
subplot(4,1,2);
plot(t,I4);
axis([0 80 0 1.08]);
ylabel('I4')
subplot(4,1,3);
plot(t,g);
axis([0 80 -1.08 1.08]);
ylabel('VT2')
subplot(4,1,4);
plot(t,h2);
axis([0 80 -1.08 1.08]);
ylabel('U2')
xlabel('Time in ms')


% VC1 = h1+f;
% VC2 = h2+g;
% 
% t = 1:length(VC1);
% tt=linspace(1,t(end),10000);
% xx=spline(t,VC1,tt);
% yy=spline(t,VC2,tt);
% plot(xx,yy,'linewidth',2)
% 
% figure
% plot(VC1,VC2,'.','markersize',12)
% ylabel('VC2')
% xlabel('VC1')
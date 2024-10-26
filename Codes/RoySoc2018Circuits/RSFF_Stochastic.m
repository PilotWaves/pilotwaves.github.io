clear

%%%External Parameters%%%
T1 = 20;  %period of R and R
T2 = 40;  %period of S and S
Total = 80;  %Total time of clock
freq = 1/10;  %internal frequency
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
Qc = f;  %U1
Q = f;  %U2
%%%%%%%%%%%%%%%%%%%%

%%%Inputs%%%
S = Qc;
S(1:T1/(2*freq)+1) = 1 - rand*1e-8;
S(T1/freq+1:3*T1/(2*freq)+1) = 1 - rand*1e-8;
S(2*T1/freq+1:5*T1/(2*freq)+1) = 1 - rand*1e-8;
S(3*T1/freq+1:7*T1/(2*freq)+1) = 1 - rand*1e-8;
R = Qc;
R(T2/(2*freq)+2:T2/freq+1) = 1 - rand*1e-8;
R(3*T2/(2*freq)+2:2*T2/freq+1) = 1 - rand*1e-8;
Qc(1) = rand*1e-8;
Q(1) = 1 - rand*1e-8;
%%%%%%%%%%%%

for n = 1:N
    
      %%%Adding noise into inputs%%%
      R(n+1) = R(n+1) + rand*1e-8;
      S(n+1) = S(n+1) + rand*1e-8;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      
      %%%VT1 calculations%%%
    
    if  f(n) >= -nu*mu && f(n) < -nu
        
        yf = (1+mu)*f(n)/(1-mu) + nu*(1+mu)/(1-mu) - mu*nu;
        
    else if f(n) >= -nu && f(n) <= nu
            
            yf = mu*f(n);
            
        else
            
            yf = (1+mu)*f(n)/(1-mu) - nu*(1+mu)/(1-mu) + mu*nu;
            
        end
        
    end
    
    f(n+1) = abs(S(n+1) - Q(n)) + S(n+1)*Q(n)*yf;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%VT2 calculations%%%
    
    if abs(g(n)-1/3) <= 1/2
    
        yg = 1/3 + 2*abs(g(n)-1/3);
        
    else
        
        yg = 1/3 + 2*(1-abs(g(n)-1/3));
        
    end
    
    g(n+1) = abs(R(n+1) - Qc(n)) - 1 + R(n+1)*Qc(n)*yg;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%RSFF Output calculations%%%
    
        Qc(n+1) = (1-S(n+1))*(1-Q(n));
        Q(n+1) = (1-R(n+1))*(1-Qc(n));
        
        T = 16;
        
     if n >= 90 && randi(2) == 1
        w1 = randperm(T,T);
        w2 = randperm(T,T);
        x1 = dot(w1,R(n-T+2:n+1))/sum(w1);
        x2 = dot(w2,S(n-T+2:n+1))/sum(w2);
%         w3 = randperm(T,T);
%         w4 = randperm(T,T);
%         x3 = dot(w3,Qc(n-T+1:n))/sum(w3);
%         x4 = dot(w4,Q(n-T+1:n))/sum(w4);
x3 = Qc(n);
x4 = Q(n);
        Qc(n+1) = (1-x2)*(1-x4);
        Q(n+1) = (1-x1)*(1-x3);
     end
    
%      if abs(R(n+1) - 1) < 1e-8 && abs(S(n+1) - 1) < 1e-8
%            Q(n+1) = g(n+1);
%            Qc(n+1) = f(n+1)-1;
%      end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

set(0,'defaultlinelinewidth',2)
set(0,'DefaultAxesFontSize',24)

subplot(4,1,1)
plot(t,S)
axis([0 80 0 1.08]);
ylabel('Set')
subplot(4,1,2)
plot(t,R)
axis([0 80 0 1.08]);
ylabel('Reset')
subplot(4,1,3)
plot(t,Q)
axis([0 80 -1.08 1.08]);
ylabel('Q')
subplot(4,1,4)
plot(t,Qc)
axis([0 80 -1.08 1.08]);
ylabel('Q^{\prime}')
xlabel('Time in ms')

figure

subplot(2,1,1)
plot(t,f);
axis([0 80 -2.16 2.16]);
ylabel('VT1')
subplot(2,1,2)
plot(t,g);
axis([0 80 -1.08 1.08]);
ylabel('VT2')
clear all

%%%External Parameters%%%
T1 = 20;  %period of I1 and I4
T2 = 40;  %period of I2 and I3
Total = 80;  %Total time of clock
freq = 1/9;  %internal frequency
%%%%%%%%%%%%%%%%%%%%%%%%%
N = Total/freq;
t = 0:freq:Total;

%%%Initialization%%%
h1 = zeros(1,N+1);  %U1
h2 = h1;  %U2
%%%%%%%%%%%%%%%%%%%%

%%%Inputs%%%
I1 = h1;
I1(1:T1/(2*freq)+1) = 1;
I1(T1/freq+1:3*T1/(2*freq)+1) = 1;
I1(2*T1/freq+1:5*T1/(2*freq)+1) = 1;
I1(3*T1/freq+1:7*T1/(2*freq)+1) = 1;
I2 = h1;
I2(T2/(2*freq)+2:T2/freq+1) = 1;
I2(3*T2/(2*freq)+2:2*T2/freq+1) = 1 + 1e-9;
I3 = I2;
I4 = I1;
h1(1) = 0;
h2(1) = 1;
%%%%%%%%%%%%

for n = 1:N
    
    %%%RSFF Output calculations%%%
    
        h1(n+1) = (1-I1(n+1))*(1-h2(n));
        h2(n+1) = (1-I3(n+1))*(1-h1(n));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

set(0,'defaultlinelinewidth',4)
set(0,'DefaultAxesFontSize',24)

subplot(4,1,1)
plot(t,I1)
axis([0 80 0 1.08]);
ylabel('Set')
subplot(4,1,2)
plot(t,I3)
axis([0 80 0 1.08]);
ylabel('Reset')
subplot(4,1,3)
plot(t,h2)
axis([0 80 -1.08 1.08]);
ylabel('Q')
subplot(4,1,4)
plot(t,h1)
axis([0 80 -1.08 1.08]);
ylabel('Q^{\prime}')
xlabel('Time in ms')
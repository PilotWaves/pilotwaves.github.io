function [Percent_kill] = EthanolDiff_CN(Vol,U_thresh,eps)

%Ethanol diffusion in tumor

%     Vol = 0.2;
%  U_thresh = 0.25;
%  eps=1;

dr = 0.01;
dt = dr^2/2;
nu = dt/dr;
mu = dt/dr^2;
T = 2;
X = 1;
N = round(T/dt + 1);
M = round(X/dr + 1);

M_kill = zeros(N,1);
Percent_kill = 0;


r = [0:dr:1];

U = zeros(M,N);
A = zeros(M);
B = A;
width = 1/4;
r1 = round(length(r)*width);
height = Vol/(sum((r(1:r1 - 1).^2).*exp(1 - 1./(1 - (r(1:r1 - 1)/r(r1)).^2)))*dr*4*pi);
U(1:r1 - 1,1) = exp(1 - 1./(1 - (r(1:r1 - 1)/r(r1)).^2))*(height);

for n = 2:N
    
    if U(1,n-1) >= U_thresh
      M_kill(n-1) = r(1)^3;
     end
    
    %Point at origin
     A(1,1) = mu + 1; A(1,2) = -mu; B(1,1) = -mu + 1; B(1,2) = mu;
     A(M,M) = mu + 1 + eps*dr*(nu + mu); A(M,M-1) = -mu;
     B(M,M) = 1 - mu - eps*dr*(nu + mu); B(M,M-1) = mu;
    
   %Points in the middle
   for m = 2:M-1
        A(m,m+1) = -nu/(2*r(m)) - mu/2; A(m,m) = mu + 1;
        A(m,m-1) = nu/(2*r(m)) - mu/2;
        B(m,m+1) = nu/(2*r(m)) + mu/2; B(m,m) = 1 - mu;
        B(m,m-1) = -nu/(2*r(m)) + mu/2;
        
           if U(m,n-1) >= U_thresh
                M_kill(n-1) = r(m)^3;
           end
   
   end
   
   U(:,n) = (B/A)*U(:,n-1);
   
      if U(M,n-1) >= U_thresh
        M_kill(n-1) = 1;
         break
      end
end

for m = 1:M
    if U(m,N) >= U_thresh
       M_kill(N) = r(m)^3;
    end
end

      Percent_kill = max(M_kill);
function [Percent_kill] = EthanolDiff_FTCS(Vol,U_thresh,eps)

%Ethanol diffusion in tumor

%     Vol = 0.2;
%  U_thresh = 0.25;

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

%NT = round(N*Time/T);
%U_leak = 0;


r = [0:dr:1];

U = zeros(N,M);
width = 1/4;
r1 = round(length(r)*width);
height = Vol/(sum((r(1:r1 - 1).^2).*exp(1 - 1./(1 - (r(1:r1 - 1)/r(r1)).^2)))*dr*4*pi);
U(1,1:r1 - 1) = exp(1 - 1./(1 - (r(1:r1 - 1)/r(r1)).^2))*(height);

for n = 2:N
    
    %Point at origin
     U(n,1) = U(n-1,1) + 2*mu*(U(n-1,2) - U(n-1,1));
     
     if U(n,1) >= U_thresh
      M_kill(n) = r(1)^3;
     end
    
   %Points in the middle
   for m = 2:M-1
      U(n,m) = U(n-1,m) + (nu/r(m))*(U(n-1,m+1) - U(n-1,m-1)) + ...
               mu*(U(n-1,m+1) - 2*U(n-1,m) + U(n-1,m-1));
%       U(n,m) = U(n-1,m) + (mu/r(m)^2)*(U(n-1,m+1)*((r(m)+r(m+1))^2)/4 ...
%                 -2*(r(m)^2)*U(n-1,m) + U(n-1,m-1)*((r(m)+r(m-1))^2)/4);
%         U(n,m) = U(n-1,m) + (mu/r(m)^2)*((U(n-1,m+1)-U(n-1,m))*((r(m)+r(m+1))^2)/4 ...
%                  - (U(n-1,m)-U(n-1,m-1))*((r(m)+r(m-1))^2)/4);

    if U(n,m) > U(n,m-1)
       blah = U(n,m-1);
       U(n,m-1) = U(n,m);
       U(n,m) = blah;
    end
            
   if U(n,m) >= U_thresh
      M_kill(n) = r(m)^3;
   end
   end
   
   %Boundary Condition
    
%     if U(n-1,M) < U_leak
%     %if n < NT+1
%        U(n,M) = U(n-1,M) + 2*mu*(U(n-1,M-1) - U(n-1,M));
%     else
       %U(n,M) = (U(n,M-2)+2*eps*dr*U_leak)/(1+2*eps*dr);
       U(n,M) = (1 - 2*dt*eps - 2*mu*(1+eps*dr))*U(n-1,M) + 2*mu*U(n-1,M-1);
       %U(n,:) = smooth(U(n,:));
       %sum((r.^2).*U(n,:))*dr
%    end
    
%         if n == NT
%             U_leak = U(NT,M);
%         end
   
   %Break Condition
   
   
   if U(n,M) >= U_thresh
        M_kill(n) = 1;
      break 
   end
      
end

      Percent_kill = max(M_kill);
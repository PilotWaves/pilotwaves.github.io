omega = 31/2;
K = -(pi*exp(omega^2/8.4))/(sin(pi*omega));
nu = omega^2/2.1/4/pi^2;

set(0,'defaultAxesFontSize',30)

% C = 1/K/10;
% C = 1/K/6;
 C = 1/K/2;

E_loss_in = 0.001*C*K;

N = 600;
E = zeros(1,N);
E(1) = E_loss_in;
E(2) = E(1)*C^2;

for i = 3:2:N-1
E(i) = E(i-1) + C^2*K^2*sin(omega*sqrt(E(i-1))/C)^2*exp(-2*nu*E(i-1)/C^2)...
             + C*K*sqrt(E(i-1))*sin(omega*sqrt(E(i-1))/C)*exp(-nu*E(i-1)/C^2);
E(i+1) = C^2*E(i);
end

nexttile
plot([501:2:N-1], E(501:2:N-1) - E(502:2:N))
if C == 1/K/10
   ytickformat('%.1f');
   ylabel('Energy Gain');
elseif C == 1/K/6
    axis([501 N min(E(501:2:N-1) - E(502:2:N)) max(E(501:2:N-1) - E(502:2:N))]);
else
    axis([501 N min(E(501:2:N-1) - E(502:2:N)) max(E(501:2:N-1) - E(502:2:N))]);
end

    xlabel('Impacts')
tic

set(0,'DefaultAxesFontSize',20)

omega = 31/2;
K = -(pi*exp(omega^2/8.4))/(sin(pi*omega));
nu = omega^2/2.1/4/pi^2;

%C = [1/K/100:1/K/10000:1/K];
C = [0.14/K:1/K/20000:1/K/5];
N = 100;
E = zeros(length(C),N);

for i = 1:length(C)
   E_loss_in = 0.1;
%    for n = 1:500
%       E_loss_out = C(i)^2*E_loss_in;
%       E_gain_in = E_loss_out;
%       E_gain_out = E_gain_in + C(i)^2*K^2*sin(omega*sqrt(E_gain_in)/C(i)).^2.*exp(-2*nu*E_gain_in/C(i)^2)...
%              + C(i)*K*sqrt(E_gain_in).*sin(omega*sqrt(E_gain_in)/C(i)).*exp(-nu*E_gain_in/C(i)^2); 
%       E_loss_in = E_gain_out;
%    end
   for n = 1:N
      E_loss_out = C(i)^2*E_loss_in;
      E_gain_in = E_loss_out;
      E_gain_out = E_gain_in + C(i)^2*K^2*sin(omega*sqrt(E_gain_in)/C(i)).^2.*exp(-2*nu*E_gain_in/C(i)^2)...
             + C(i)*K*sqrt(E_gain_in).*sin(omega*sqrt(E_gain_in)/C(i)).*exp(-nu*E_gain_in/C(i)^2); 
      E_loss_in = E_gain_out;
      E(i,n) = E_gain_out;
   end
      if round(E(i,N),4) == round(E(i,N-1),4)
          Col = [0.1, 0.1, 0.1];
          plot(C(i)*K,E(i,:),'.','Color',Col, 'MarkerSize', 1)
          hold on
      elseif round(E(i,N),4) == round(E(i,N-2),4)
          Col = [0, 0.4470, 0.7410];
          plot(C(i)*K,E(i,:),'.','Color',Col, 'MarkerSize', 2)
          hold on
      elseif round(E(i,N),4) == round(E(i,N-4),4)
          Col = [0.9290, 0.6940, 0.1250];
          plot(C(i)*K,E(i,:),'.','Color',Col, 'MarkerSize', 4)
          hold on
      elseif round(E(i,N),4) == round(E(i,N-8),4)
          Col = [0.4940 0.1840 0.5560];
          plot(C(i)*K,E(i,:),'.','Color',Col, 'MarkerSize', 6)
          hold on
      elseif round(E(i,N),4) == round(E(i,N-3),4)
          Col = [1, 0, 0];
          plot(C(i)*K,E(i,:),'.','Color',Col, 'MarkerSize', 10)
          hold on
      else
          %if randi([0, 1]) == 1
          Col = [0.4660, 0.6740, 0.1880];
          %else
          %    Col = [0.3010, 0.7450, 0.9330];
          %end
          plot(C(i)*K,E(i,:),'.','Color',Col, 'MarkerSize', 0.5)
          hold on
      end
end

hold off
xticks([0.15, 0.175, 0.2])
xticklabels({'3/20K', '7/40K', '1/5K'})
ylabel('E_{total}')
xlabel('Bifurcation Parameter C')


%plot(C*K,E,'.', 'Color', [0.2, 0.2, 0.2], 'MarkerSize', 0.5)
%xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
%xticklabels({'0', '1/5K', '2/5K', '3/5K', '4/5K', '1/K'})
%ylabel('E_{total}')
%xlabel('Bifurcation Parameter C')
toc
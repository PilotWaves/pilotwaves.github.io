omega = 31/2;
K = -(pi*exp(omega^2/8.4))/(sin(pi*omega));
C = 1/K;
nu = omega^2/2.1/4/pi^2;


E_loss_in = [0:0.001*C*K:pi^2];
E_gain_in = E_loss_in*C^2;

E_loss_out = C^2*E_loss_in;
E_gain_out = E_gain_in + C^2*K^2*sin(omega*sqrt(E_gain_in)/C).^2.*exp(-2*nu*E_gain_in/C^2)...
             + C*K*sqrt(E_gain_in).*sin(omega*sqrt(E_gain_in)/C).*exp(-nu*E_gain_in/C^2);

E_loss_in = [min(E_gain_out):0.001*C*K:max(E_gain_out)*1.1];
E_loss_out1 = C^2*E_loss_in;
E_loss_out4 = (C/4)^2*E_loss_in;
E_loss_out20 = (C/20)^2*E_loss_in;

set(0,'defaultAxesFontSize',20)
t = tiledlayout(1,2,'TileSpacing','Compact');

plot(E_loss_in, E_loss_out1, '--', E_gain_out, E_gain_in, E_loss_in, E_loss_out4, '--', E_loss_in, E_loss_out20, '--', 'linewidth', 4)
axis([min(E_loss_in) max(E_loss_in) min(E_loss_out1) max(E_loss_out1)/3])
xlabel({'E_{in}^{loss}, E_{out}^{gain}'})
ylabel({'E_{out}^{loss}, E_{in}^{gain}'})

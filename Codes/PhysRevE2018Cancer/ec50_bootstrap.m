 clear all
 r=[0 0.0406 0.1075 0.4026 0.7591 0.8621]; % response
 d=[0 0.0316 0.1 0.316 1 3.16];% dose
%function ec50_bootstrap(d,r)

       [coeff sim_resid]=sigmoid_doseResponse(d,r);
       A_min=coeff(1);
       A_max=coeff(2);
       ec50=coeff(3);
       hill=coeff(4);
       
       mean_fn1=@(coeff,d) coeff(1)+(coeff(2)-coeff(1))./(1+(d/coeff(3)).^coeff(4));
       mean_fn=[mean_fn1(coeff,d)]';
       %plot(d,r,'o'); hold on, plot(d, mean_fn)

       nboot=50;
       sigma= 0.0003;theta=12; %% can tune these 2 parameters to reduce noise
       [bb1,bb2]=meshgrid(d);
       Sigma=sigma*exp(-theta*(bb2-bb1).^2);
       rng(20)
       gp_err=[mvnrnd(zeros(1,length(d)),Sigma,nboot)]';
       yy_boot=repmat(mean_fn,1,nboot)+ gp_err;
       %plot(d,yy_boot)
       
%        ec50_boot=[];
%        for ii=1:nboot
%            [coeff_boot sim_resid_boot]=sigmoid_doseResponse(d',yy_boot(:,ii));
%            ec50_boot=[ec50_boot;(coeff_boot(3))];
%        end
%quant = quantile(ec50_boot,[.025,.5,.975]);
quant = quantile(yy_boot,[.025,.975],2);

set(0,'defaultAxesFontSize',14)
%plot(d,quant(:,1),'s',d,quant(:,2),'s','MarkerSize',20)
plot(d,quant(:,1),'--',d,quant(:,2),'--','LineWidth',2)
%hold on

% plot(d,yy_boot)
% xlabel('Initial Average Concentration (Fraction of Tumor Volume)')
% ylabel('Percent Cells Killed')
% alw = 0.75;    % AxesLineWidth
% axis([0 max(d) 0 max(r)])

% plot([min(d) max(d)],[max(r)/2 max(r)/2],'--','LineWidth',2)
% axis([min(d) max(d) min(r) max(r)])
% xlabel('Initial Average Concentration (Fraction of Tumor Volume)')
% ylabel('Apoptosis fraction')
% alw = 0.75;    % AxesLineWidth
% title({'Cell line: C32','Drug: Selumetinib'})
% hold off
% txt = ['[',num2str(quant(1)),',',num2str(quant(3)),']'];
% annotation('textbox',[0.5 0.3 0.9 0.1],'FontSize',14,'String',txt,'FitBoxToText','on');
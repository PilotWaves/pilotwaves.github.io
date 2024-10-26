N = 100000;
t = [1:N];

v = zeros(1,N);
theta = v;

v(1) = 0.1;
theta(2) = 1;

set(0,'DefaultAxesFontSize',14)

K = 8;

for n = 2:N-1
    
    v(n+1) = mod(v(n) + K*exp(-abs(theta(n)-theta(n-1)))*sin(theta(n)),2*pi);
    % v(n+1) = mod(v(n) + K*sin(theta(n)),2*pi);
    theta(n+1) = mod(theta(n) + v(n+1),2*pi);
    
end


% h = plot(t,theta,'.','MarkerSize',20);
% xlabel('n')
% ylabel('\theta')

h = plot(theta,v,'.','MarkerSize',1);
xlabel('\theta')
ylabel('v')
axis([0 2*pi 0 2*pi])

%%% Histograms %%%
% 
% [f,theta] = hist(theta,N);
% 
% bar(theta,f/sum(f));
% alw = 0.75;    % AxesLineWidth
% fsz = 28;      % Fontsize
% xlabel('\theta')
% ylabel('Probability')

% x = cos(theta);
% y = sin(theta);
% z = f;
% col = f;  % This is the color, vary with f in this case.
% surface([x;x],[y;y],[z;z],[col;col],...
%         'facecol','no',...
%         'edgecol','interp',...
%         'linew',5);
% axis([-1.5 1.5 -1.5 1.5])
% axis off

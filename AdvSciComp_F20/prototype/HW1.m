%% Homework 1 -- Alanna Gary
close all
clear all

%% Problem 1(a): FE
f = @(t,y) -3*y*sin(t); %eqn for dy/dt
true_soln = @(t) pi.*exp(3.*(cos(t)-1))./sqrt(2); %true soln to ODE
dt_vals = 2.^(-2:-1:-8); %the values for delta t
Tmax = 5;
y_0 = pi/sqrt(2);
Err = zeros(1,length(dt_vals));

for i=1:length(dt_vals);
    t = 0:dt_vals(i):Tmax; %set t for this value of dt
    Y = zeros(1,length(t));
    Y(1) = y_0;
    for j=1:(length(t)-1)
        Y(j+1) = Y(j) + dt_vals(i)*f(t(j),Y(j)); %Do the Euler!
    end
    Err(i) = mean(abs(true_soln(t) - Y)); %Calculate 2-norm of err
%     figure;
%     plot(t,Y,t,true_soln(t),'--');
%     legend('Approx','True');
%     xlabel('Time');ylabel('Soln');
%     title(['FE solns, \Delta t = ',num2str(dt_vals(i))]);
end
ans_1 = transpose(Y); %make column vector
save('A1.dat','ans_1','-ascii');
save('A2.dat','Err','-ascii'); %save err as row vector
 
% figure;
% plot(log(dt_vals),log(Err));
% xlabel('log(\Delta t)');ylabel('log(Error)');
% title('FE error vs. \Delta t log-log plot');
p = polyfit(log(dt_vals),log(Err),1);
slope = p(1);
save('A3.dat','slope','-ascii');



%% Problem 1(b): Heun's method

Err = zeros(1,length(dt_vals));
for i=1:length(dt_vals);
    t = 0:dt_vals(i):Tmax; %set t for this value of dt
    Y = zeros(1,length(t));
    Y(1) = y_0;
    for j=1:(length(t)-1) %Do the Heun!
        Y(j+1) = Y(j) + dt_vals(i)/2*(f(t(j),Y(j)) + f(t(j+1), Y(j) + dt_vals(i)*f(t(j),Y(j)))); %Do the Euler!
    end
    Err(i) = mean(abs(true_soln(t) - Y)); %Calculate 2-norm of err
%     figure;
%     plot(t,Y,t,true_soln(t),'--');
%     legend('Approx','True');
%     xlabel('Time');ylabel('Soln');
%     title(['Heun solns, \Delta t = ',num2str(dt_vals(i))]);
end
ans_4 = transpose(Y); %make column vector
save('A4.dat','ans_4','-ascii');
save('A5.dat','Err','-ascii');

% figure;
% plot(log(dt_vals),log(Err));
% xlabel('log(\Delta t)');ylabel('log(Error)');
% title('Heun error vs. \Delta t log-log plot');
p = polyfit(log(dt_vals),log(Err),1);
slope = p(1);
save('A6.dat','slope','-ascii');



%% Problem 2(a): VDP ode45
%First, break into system of 2 first-order ODEs
%Let x1 = y, x2 = y'. Then x1' = x2, and x2' = -eps*(x1^2 - 1]*x2 - x1
f = @(t,x,eps) [x(2); -eps*(x(1)^2 - 1)*x(2) - x(1)];
x_init = [sqrt(3); 1];
eps_vals = [0.1 1 20]; %set values of epsilon
t = 0:0.5:32;
Ans = zeros(length(t),length(eps_vals));
for i=1:length(eps_vals);
    [t,X] = ode45(@(t,x) f(t,x,eps_vals(i)),t,x_init);
    Ans(:,i) = X(:,1); %y is saved as x(:,1); y' = x(:,2).
end
save('A7.dat','Ans','-ascii');



%% Problem 2(b): VDP comparison of methods
x_init = [2; pi^2];
eps = 1;
t = [0 32];
tol_vals = 10.^(-4:-1:-10);
o45_meantime = zeros(length(tol_vals),1);
o23_meantime = zeros(length(tol_vals),1);
o113_meantime = zeros(length(tol_vals),1);

for i=1:length(tol_vals)
    tol = tol_vals(i);
    options = odeset('RelTol',tol,'AbsTol',tol);
    %Ode45
    [time,~] = ode45(@(t,x) f(t,x,eps),t,x_init,options);
    o45_meantime(i) = mean(diff(time));
    
    %Ode23
    [time,~] = ode23(@(t,x) f(t,x,eps),t,x_init,options);
    o23_meantime(i) = mean(diff(time));
    
    %Ode113
    [time,~] = ode113(@(t,x) f(t,x,eps),t,x_init,options);
    o113_meantime(i) = mean(diff(time));
end


% figure;
% plot(log(o45_meantime),log(tol_vals));
% xlabel('log(tol)');ylabel('log(\Delta t)');
% title('LTE for ODE45');
p = polyfit(log(o45_meantime),log(tol_vals)',1);
slope45 = p(1);

% figure
% plot(log(o23_meantime),log(tol_vals));
% xlabel('log(tol)');ylabel('log(\Delta t)');
% title('LTE for ODE23');
p = polyfit(log(o23_meantime),transpose(log(tol_vals)),1);
slope23 = p(1);

% figure
% plot(log(o113_meantime),log(tol_vals));
% xlabel('log(tol)');ylabel('log(\Delta t)');
% title('LTE for ODE113');
p = polyfit(log(o113_meantime),transpose(log(tol_vals)),1);
slope113 = p(1);

save('A8.dat','slope45','-ascii');
save('A9.dat','slope23','-ascii');
save('A10.dat','slope113','-ascii');

%% Problem 3: coupled FH model
a1 = 0.05;
a2 = 0.25;
b = 0.01;
c = 0.01;
I = 0.1;
d12 = -0.1; %NEGATIVE (arbitrary)
d21 = 0.1; %POSITIVE (arbitrary)

param_vec(1) = a1;
param_vec(2) = a2;
param_vec(3) = b;
param_vec(4) = c;
param_vec(5) = I;
param_vec(6) = d12;
param_vec(7) = d21;

init_vec = [0.1;0;0.1;0]; %set initial condition, [v1;w1;v2;w2]

%Following is asked for, but not really graded?:
% [t,Y] = ode15s(@(t,in_vec) coupledFH(t,in_vec,param_vec), 0:0.5:100, init_vec);
% figure;
% plot(Y(:,1),Y(:,2));
% xlabel('v1');ylabel('w1');
% title('Phase portrait, first neuron');
% figure;
% plot(Y(:,3),Y(:,4));
% xlabel('v2');ylabel('w2');
% title('Phase portrait, second neuron');
% figure;
% plot(t,Y(:,1),t,Y(:,3));
% xlabel('Time');ylabel('Voltage');
% title('Time-plot of voltage, both neurons');
% legend('Neuron 1','Neuron 2');
% With these interaction parameters, the neurons appear to fire slightly
% out-of-phase with one another.


%Keeping init_vec the same; changing d12, d21 values.
d_vec = [0 0; 0 0.2; -0.1 0.2; -0.3 0.2; -0.5 0.2]; %(d12,d21)
t = 0:0.5:100;
big_answer_array = zeros(length(t),4,length(d_vec));
for i=1:length(d_vec)
    param_vec(6) = d_vec(i,1); %d12
    param_vec(7) = d_vec(i,2); %d21
    [~,Y] = ode15s(@(t,in_vec) coupledFH(t,in_vec,param_vec), t, init_vec);
    big_answer_array(:,1,i) = Y(:,1); %v1
    big_answer_array(:,2,i) = Y(:,3); %v2
    big_answer_array(:,3,i) = Y(:,2); %w1
    big_answer_array(:,4,i) = Y(:,4); %w2
end

soln1 = big_answer_array(:,:,1);
save('A11.dat','soln1','-ascii');

soln2 = big_answer_array(:,:,2);
save('A12.dat','soln2','-ascii');

soln3 = big_answer_array(:,:,3);
save('A13.dat','soln3','-ascii');

soln4 = big_answer_array(:,:,4);
save('A14.dat','soln4','-ascii'); 

soln5 = big_answer_array(:,:,5);
save('A15.dat','soln5','-ascii');

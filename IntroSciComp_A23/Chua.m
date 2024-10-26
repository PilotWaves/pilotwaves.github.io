[T,Y] = ode45(@(t,x) ChuaODE(t,x), [0:0.05:100], [0.1; 0.2; 0.3]);
plot3(Y(:,1), Y(:,2), Y(:,3))
%plot(T, Y(:,1))

function dx = ChuaODE(t,x)

alpha = 16;
beta = 100;
    dx1 = alpha*(x(2) - (x(1)^3/16 - x(1)/6));
    dx2 = x(1) - x(2) + x(3);
    dx3 = -beta*x(2);
    dx = [dx1; dx2; dx3];
end
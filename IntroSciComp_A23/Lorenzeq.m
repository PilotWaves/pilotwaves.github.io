function dx = Lorenzeq(t, x, params)

sigma = params(1);
beta = params(2);
rho = params(3);

dx1 = -sigma*x(1) + sigma*x(2);
dx2 = rho*x(1) - x(2) - x(1)*x(3);
dx3 = -beta*x(3) + x(1)*x(2);

dx = [dx1; dx2; dx3];

end
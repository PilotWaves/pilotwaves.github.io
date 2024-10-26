%InjIC1:  The IC is given as a decaying exponential that jumps to zero
%         after injection stops.

alpha = 4; %  internal diffusion constant

r = [0:.01:1];
t = [0:.01:4];

Psi0 = 1;

A = Psi0;
u = exp(-t)*(A*cos(r/sqrt(alpha)) + B*sin(r/sqrt(alpha)));

x = r'*cos([0:pi/50:2*pi]);
y = r'*sin([0:pi/50:2*pi]);
z = u(1,:)'*ones(size(r));
% surf(x,y,u)
col = z;  % This is the color, vary with u in this case.
surface(x,y,z,col,...
        'facecol','no',...
        'edgecol','interp',...
        'linew',5);
axis([-1.5 1.5 -1.5 1.5])
axis off
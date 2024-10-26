function dy = DoublePendulum(t,y)
dy = zeros(4,1);    % a column vector

m1 = 2; %mass of the first in kilograms
m2 = 1; %mass of the second in kilograms
L1 = 2; %length of the first in meters
L2 = 1; %length of the second in meters
g = 9.8; %gravity

%y(1): phi_1
%y(2): phi_2
%y(3): phi_1'
%y(4): phi_2'

dy(1) = y(3); %phi_1'

dy(2) = y(4); %phi_2'

dy(3) = (1./((m1+m2)*L2 - m2*L1*(cos(y(1) - y(2))).^2)).*(-m2*L2*cos(y(1) - y(2)).*((y(3) - y(4)).*L1.*y(3).*sin(y(1) - y(2)) + L1*y(3).*y(4).*sin(y(1) - y(2)) + g*sin(y(2))) + (y(3) - y(4)).*m2.*L2.*y(4).*sin(y(1) - y(2)) - m2*L2*y(3).*y(4).*sin(y(1) - y(2)) - (m1 + m2)*g.*sin(y(1))); %theta_1'                                                  

dy(4) = (L1/L2)*(((y(3) - y(4)).*L1.*y(3).*sin(y(1) - y(2)) + L1*y(3).*y(4).*sin(y(1) - y(2)) + g*sin(y(2))) - (cos(y(1) - y(2))./((m1+m2)*L2 - m2*L1*(cos(y(1) - y(2))).^2)).*(-m2*L2*cos(y(1) - y(2)).*((y(3) - y(4)).*L1.*y(3).*sin(y(1) - y(2)) + L1*y(3).*y(4).*sin(y(1) - y(2)) + g*sin(y(2))) + (y(3) - y(4)).*m2.*L2.*y(4).*sin(y(1) - y(2)) - m2*L2*y(3).*y(4).*sin(y(1) - y(2)) - (m1 + m2)*g.*sin(y(1)))); %theta_2'
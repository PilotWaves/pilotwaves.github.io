function dy = SimplePendulum(t,y)
    dy = zeros(2,1);
    
    dy(1) = y(2);
    dy(2) = -sin(y(1));
end
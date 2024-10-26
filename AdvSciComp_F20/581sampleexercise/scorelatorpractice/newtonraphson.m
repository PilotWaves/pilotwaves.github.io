% Newton-Raphson algorithm to find a root

function [xvals,steps]=newtonraphson(x0,tol)

%define the input parameters


xn = x0;

steps = 0;
xvals(1,1) = xn;

N = 1000;

% iterate NR scheme

for i=1:N
    
   % calculate f and its derivative 
   
   steps = steps+1;
   
   fxn = func(xn);
   fprimexn = funcder(xn); 
   xnplus1 = xn - fxn/fprimexn;
   
   fxnplus1 = func(xnplus1);
   
   xvals(end+1,1) = xnplus1;
      
   if (abs(fxnplus1) < tol)

        break;
   end
   
   xn = xnplus1;

end
% Bisection algorithm to find a root

function [midvals,steps]=bisection(a,b,tol)

%define the input parameters


N = 1000;

% for testing purposes

steps = 0;
midvals= [];
for j=1:N
    
   % increse steps 
   steps = steps+1;
    
   % chose a midpoint 
   c = (a+b)/2;
   
   % save the midpoint value in the last row of midvals
   midvals(end+1,1) = c;
   
   % evaluate f at c
   fc = func(c);
   
   % plot the guess
   %hold on;
   %plot(c,fc,'*r')
   
   % adjust the interval according to f(c)
   if (fc > 0)
       a=c;
       
   else
       b=c;
   end
   
    
   if (abs(fc) < tol)

        break;
   end
    
end

% Homework 2 MATLAB template file
% Your main file should be named "solution.m" and it should be saved as UTF-8 file.

function [consoleout, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13] = solution()
 [consoleout, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13] = evalc('student_solution(0)'); 
end

function [A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13] = student_solution(dummy_argument)
    % your solution code goes here 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Problem 1
    A1 = [34, 45; 17, 6]; %Solution to Problem 1
    
    %Problem 2
    A = [1, 2; -1, 1];
    B = [2, 0; 0, 2];
    C = [2, 0, -3; 0, 0, -1];
    D = [1, 2; 2, 3; -1, 0];
    x = [1; 0];
    y = [0; 1];
    z = [1; 2; -1];
    
    A2 = A+B; %Solution to Problem 2a
    A3 = 3*x - 4*y; %Solution to Problem 2b
    A4 = A*x; % Solution to Problem 2c
    A5 = B*(x-y); %Solution to Problem 2d
    A6 = D*x; %Solution to Problem 2e
    A7 = D*y + z; %Solution to Problem 2f
    A8 = A*B; %Solution to Problem 2g
    A9 = B*C; %Solution to Problem 2h
    A10 = C*D; %Solution to Problem 2i
    
    %Problem 3
    [A11, A13_1] = Newton_Ralphson(-3);
    
    [A12, A13_2] = Bisection(-3, 1);
    
    A13 = [A13_1, A13_2];
    
end

% your extra functions, if you need them, can be here or in another file (don't forget to upload them too!)
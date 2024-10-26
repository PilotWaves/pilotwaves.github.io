function A = GE(A)
    n = size(A); % gives a vector [n m] when A is n x m
    n = n(1); % n is now the number of rows
    for i = 1:n-1
       if A(i,i) == 0.  % Check if A(i,i) is zero
          A = 'Method failed: matrix is rank deficient';
          break;
       end
       for j = i+1:n % perform the row reduction
          A(j,:) = A(j,:) - A(j,i)/A(i,i)*A(i,:); 
       end    
    end
    
end
% Gauss Seidel method, AX=B

A = input('Enter the coefficient matrix A: ');
B = input('Enter the constant matrix B: ');

P=[A B]; % augmented matrix
[row, col] = size(P);

X=zeros(row,1);
temp=zeros(row,1);
Err=ones(row,1);

for m = 1:row % if strictly diagonally dominant matrix
    if 2 * abs(A(m,m)) <= sum(abs(A(m,:))) 
        disp("Gauss Seidel method can't be applied.");
        return
    end
end

max_err = max(Err);
while max_err > 0.00001
    for m = 1:1:row
       temp(m,1) = X(m,1);
       X(m,1) = (1 / P(m,m)) * (P(m,col) - sum(A(m,:) * X(:,1)) + A(m,m) * X(m,1));
       Err(m,1) = abs(temp(m,1) - X(m,1));
       temp(m,1) = X(m,1);
    end
    max_err=max(Err);
end

disp(' The answer is: X = ');
disp(X);


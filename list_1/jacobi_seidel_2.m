function [x_jacobi, x_seidel] = jacobi_seidel_2(A, x_star, eps)
    b = A*x_star;
    
    % Jacobi
    [m, ~] = size(A);
    x = zeros(m, 1);
    C = A;
    C(logical(eye(size(C)))) = 0;
    C = -C;
    A_diag = diag(A);
    C = C./A_diag;
    d = b./A_diag;
    err = 1;
    while err > eps
        x = C*x + d;
        err = norm(x - x_star, 'inf');
    end
    x_jacobi = x;
    
    % Seidel
    x = zeros(m, 1);
    err = 1;
    while err > eps
        for i=1:m
            x(i) = A(i, 1:(i-1))*x(1:(i-1)) + A(i, (i+1):m)*x((i+1):m);
            x(i) = -x(i);
            x(i) = x(i) + b(i);
            x(i) = x(i)/A(i, i);
        end
        err = norm(x - x_star, 'inf');
    end
    x_seidel = x;
end

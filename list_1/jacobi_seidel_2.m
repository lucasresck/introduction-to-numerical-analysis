function [x_jacobi, x_seidel] = jacobi_seidel_2(A, b, eps)
    % Jacobi
    [m, ~] = size(A);
    x = zeros(m, 1);
    C = A;
    C(logical(eye(size(C)))) = 0;
    C = -C;
    A_diag = diag(A);
    C = C./A_diag;
    d = b./A_diag;
    C_norm = norm(C, 'inf');
    factor = C_norm / (1 - C_norm);
    x_old = x;
    err = 1;
    while err > eps
        x = C*x + d;
        err = factor*norm(x - x_old, 'inf');
        x_old = x;
    end
    x_jacobi = x;
    
    % Seidel
    alpha = 0;
    for i=1:m
        p_i = sum(abs(A(i, 1:(i-1))))/abs(A(i, i));
        q_i = sum(abs(A(i, (i+1):m)))/abs(A(i, i));
        frac = q_i/(1 - p_i);
        if frac > alpha
            alpha = frac;
        end
    end
    frac = alpha/(1 - alpha);
    x = zeros(m, 1);
    x_old = x;
    err = 1;
    while err > eps
        for i=1:m
            x(i) = A(i, 1:(i-1))*x(1:(i-1)) + A(i, (i+1):m)*x((i+1):m);
            x(i) = -x(i);
            x(i) = x(i) + b(i);
            x(i) = x(i)/A(i, i);
        end
        err = frac*norm(x-x_old, 'inf');
        x_old = x;
    end
    x_seidel = x;
end


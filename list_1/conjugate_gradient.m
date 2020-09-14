function [err_jacobi, err_seidel, x] = conjugate_gradient(n, k_jacobi, k_seidel)
    % Check if n is OK
    if n < 4
        error('n is less than 4.');
    else
        if mod(n, 2) ~= 0
            error('n is not an ever number.');
        end
    end
    % Start x(0)
    x = zeros(n, 1);
    % Start r_0 = b - Ax(0) = b = (5/2, ...)
    r = zeros(n, 1);
    r(1:n) = 3/2;
    r(1) = 5/2;
    r(n) = 5/2;
    r(n/2) = 1;
    r(n/2+1) = 1;
    % Start d_0 = r_0
    d = r;
    
    err_jacobi = 0;
    err_seidel = 0;
    
    for i=1:n
        % Stop condition
        if sum(abs(r)) == 0
            break
        end
        % Calculate the denominator of various fractions
        % d_k^T*A*d_k
        mult = zeros(n, 1);
        for j=1:n
            if j == 1
                mult(1) = 3*d(1)- d(2) +1/2*d(n);
            else
                if j == n
                    mult(n) = 3*d(n)- d(n-1) +1/2*d(1);                                        
                else
                    if j == n/2 || j == n/2 + 1
                        mult(j) = 3*d(j) - d(j-1) - d(j+1);
                    else
                        mult(j) = 3*d(j) - d(j-1) - d(j+1) + 1/2*d(n-j+1);
                    end
                end
            end
        end
        denom = d'*mult;
        % Calculate alpha_k
        alpha = d'*r/denom;
        x = x + alpha*d;
        r = r - alpha*mult;
        % Calculate the beta coefficient using old calculated numbers and
        % vectors
        beta = r'*mult/denom;
        d = r - beta*d;
        if i == k_jacobi
            err_jacobi = norm(x-1, 'inf');
        end
        if i == k_seidel
            err_seidel= norm(x-1, 'inf');
        end
        if i >= k_jacobi && i >= k_seidel
            break
        end
    end
end

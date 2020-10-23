function x_t = anchored_cubic_spline(x, k_0, k_1)
    x_hat_0 = x(2);
    x_hat_n = x(end);
    x(end) = [];
    x(2) = [];
    [n, ~] = size(x);
    n = n-1;
    h = 1/n;
    A = zeros(n+1, n+1);
    b = zeros(n+1, 1);
    for i=0:n
        if i == 0
            j = i+1;
            A(j, j) = 1/(3*n);
            A(j, j+1) = 1/(6*n);
            b(j) = -k_0*(x_hat_0 - x(j)) + n*(x(j+1) - x(j));
        else
            if i == n
                j = i+1;
                A(j, j-1) = 1/(6*n);
                A(j, j) = 1/(3*n);
                b(j) = k_1*(x_hat_n - x(j)) - n*(x(j) - x(j-1));
            else
                j = i+1;
                A(j, j-1) = 1/(6*n);
                A(j, j) = 2/(3*n);
                A(j, j+1) = 1/(6*n);
                b(j) = n*(x(j+1) - 2*x(j) + x(j-1));
            end
        end
        
    end
    % Atualizar com algoritmo próprio
    M = linsolve(A, b);
%     M_0 = -3*n*k_0*(x_hat_0-x(1)) + 3*n^2*(x(2)-x(1)) - M(1)/2;
%     M_n = 3*n*k_1*(x_hat_n-x(n+1)) - 3*n^2*(x(n+1)-x(n)) - M(n-1)/2;
%     M = [M_0; M; M_n];
    x_t = [];
    for i=0:(n-1)
        p = interpolator(M(i+1), M(i+2), i/n, (i+1)/n, x(i+1), x(i+2), n);
%         p = diff(diff(p));
        x_t = [x_t, p(i/n:0.01:(i+1)/n)];
    end
end
function interactive_spline(n, k, k_0, k_1)
    clf;
    axis([0, 10, 0, 10]);
    hold on;
    [x_click, y_click] = ginput(1);
    x = [x_click];
    y = [y_click];
    plot(x, y, 'o', ...
        'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'k');
    for j = 1:k
        for i = 1:n+2
            [x_click, y_click] = ginput(1);
            x = [x; x_click];
            y = [y; y_click];
            if i == 1 || i == n+2
                plot(x_click, y_click, 'o', ...
                'MarkerEdgeColor', 'r', ...
                'MarkerFaceColor', 'r');
                plot(x(end-1:end), y(end-1:end), '--', ...
                    'Color', [200, 200, 200]/255);
            else
                if i == n+1
                    plot(x_click, y_click, 'o', ...
                    'MarkerEdgeColor', 'k', ...
                    'MarkerFaceColor', 'k');
                else
                    plot(x_click, y_click, 'o', ...
                    'MarkerEdgeColor', 'b', ...
                    'MarkerFaceColor', 'b');
                end
            end
        end
        parametric_cubic_spline(x, y, k_0, k_1);
        x = [x(end-1)];
        y = [y(end-1)];
    end
    hold off;
end

function parametric_cubic_spline(x, y, k_0, k_1)
    x_t = anchored_cubic_spline(x, k_0, k_1);
    y_t = anchored_cubic_spline(y, k_0, k_1);
    plot(x_t, y_t, '-k');
end

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
    M = seidel(A, b, 10^-6);
    x_t = [];
    for i=0:(n-1)
        p = interpolator(M(i+1), M(i+2), i/n, (i+1)/n, x(i+1), x(i+2), n);
        x_t = [x_t, p(i/n:0.01:(i+1)/n)];
    end
end

function x = seidel(A, b, eps)
    [m, ~] = size(A);
    x = zeros(m, 1);
    err = 1;
    x_old = x;
    while err > eps
        for i=1:m
            x(i) = A(i, 1:(i-1))*x(1:(i-1)) + A(i, (i+1):m)*x((i+1):m);
            x(i) = -x(i);
            x(i) = x(i) + b(i);
            x(i) = x(i)/A(i, i);
        end
        err = norm(x - x_old, 'inf')/norm(x_old, 'inf');
        x_old = x;
    end
end

function x_t = interpolator(M_i, M_ip, t_i, t_ip, x_i, x_ip, n)
    c_1 = (x_ip-x_i)*n - (M_ip - M_i)/(6*n);
    c_2 = x_i - M_i/(6*n^2) - c_1*t_i;
    syms x_t(t);
    x_t(t) = n*M_i/6*(t_ip-t)^3 + n*M_ip/6*(t-t_i)^3 + c_1*t + c_2;
end
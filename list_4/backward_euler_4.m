function [t, x] = backward_euler_4(h, t_0, T, x_0)
    % Approximate the solution with Backward Euler method.
    % Examples:
        % [t, x] = backward_euler_4(0.0001, 0, 0.5, [10;10]);
        % [t, x] = backward_euler_4(0.01, 0, 0.5, [10;10]);
    x = [x_0];
    t = [t_0];
    N = floor((T - t_0)/h);
    for k = 1:N
        x_k_old = x(:, end);
        t_k_old = t(end);
        x_k = seidel(x_k_old, 10^-3, h);
        x = [x, x_k];
        t_k = t_0 + k*h;
        t = [t, t_k];
    end
    plot_solutions(t, x, x_0)
end

function x = seidel(b, eps, h)
    x = zeros(2, 1);
    err = 1 + eps;
    x_old = x;
    C = [0, h/(1+1000*h); 0, 0];
    d = b ./ [1+1000*h;1+1/10*h];
    while err > eps
        x = C*x_old + d;
        err = norm(x - x_old, 'inf')/norm(x_old, 'inf');
        x_old = x;
    end
end

function plot_solutions(t, x, x_0)
    % Plot the exact and numerical solutions.
    y1 = x(1, :);
    y2 = x(2, :);
    y = exact_x(t, x_0);
    tiledlayout(2,1)
    
    % Plot the first coordinate
    nexttile
    plot(t, y1)
    hold on
    plot(t, y(1, :))
    hold off
    error_1 = (norm(y1 - y(1, :)));
    title(sprintf('x_1(t), with error %0.3e', error_1))

    % Plot the second coordinate
    nexttile
    plot(t,y2)
    hold on
    plot(t, y(2, :))
    hold off
    error_2 = (norm(y2 - y(2, :)));
    title(sprintf('x_2(t), with error %0.3e', error_2))
end

function x = exact_x(t, x_0)
    % Calculate the exact value of x(t).
    c_1 = x_0(1);
    c_2 = x_0(2);
    x_1 = c_1*exp(-1000*t) + c_2*(exp(-t/10)/1000 - exp(-1000*t)/1000);
    x_2 = c_2*exp(-t/10);
    x = [x_1; x_2];
end


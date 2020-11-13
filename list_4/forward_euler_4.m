function [t, x] = forward_euler_4(h, t_0, T, x_0)
    % Approximate the solution with Forward Euler method.
    % Examples:
        % [t, x] = forward_euler_4(0.0001, 0, 0.5, [10;10]);
        % [t, x] = forward_euler_4(0.01, 0, 0.5, [10;10]);
    x = [x_0];
    t = [t_0];
    N = floor((T - t_0)/h);
    for k = 1:N
        x_k_old = x(:, end);
        t_k_old = t(end);
        x_k = x_k_old + f(x_k_old)*h;
        x = [x, x_k];
        t_k = t_0 + k*h;
        t = [t, t_k];
    end
    plot_solutions(t, x, x_0, h)
end

function result = f(x)
    % Calculate f.
    A = [-1000, 1; 0, -1/10];
    result = A*x;
end

function plot_solutions(t, x, x_0, h)
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
    legend({'Numerical solution', 'Exact solution'}, ...
    'Location', 'northeast')
    error_1 = (norm(y1 - y(1, :)));
    title(sprintf('x_1(t), h = %0.5f, error = %0.3e', h, error_1))

    % Plot the second coordinate
    nexttile
    plot(t,y2)
    hold on
    plot(t, y(2, :))
    hold off
    legend({'Numerical solution', 'Exact solution'}, ...
    'Location', 'northeast')
    error_2 = (norm(y2 - y(2, :)));
    title(sprintf('x_2(t), h = %0.5f, error = %0.3e', h, error_2))
end

function x = exact_x(t, x_0)
    % Calculate the exact value of x(t).
    c_1 = x_0(1);
    c_2 = x_0(2);
    x_1 = c_1*exp(-1000*t) + c_2*(exp(-t/10)/1000 - exp(-1000*t)/1000);
    x_2 = c_2*exp(-t/10);
    x = [x_1; x_2];
end

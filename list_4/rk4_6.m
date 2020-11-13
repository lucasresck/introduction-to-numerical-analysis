function x = rk4_6(h, T, k, a, b)
    % Calculate the Runge–Kutta fourth-order method
    % Examples:
        % x = rk4_6(0.5, 20, 0.01, 70, 50);
    x_0 = 0;
    t_0 = 0;
    x = [x_0];
    t = [t_0];
    N = floor((T - 0)/h);
    for i = 1:N
        x_i_old = x(end);
        k_1 = f(x_i_old, k, a, b);
        k_2 = f(x_i_old + h/2*k_1, k, a, b);
        k_3 = f(x_i_old + h/2*k_2, k, a, b);
        k_4 = f(x_i_old + h*k_3, k, a, b);
        phi = 1/6*(k_1 + 2*k_2 + 2*k_3 + k_4);
        x_i = x_i_old + h*phi;
        x = [x, x_i];
        t_i = t_0 + i*h;
        t = [t, t_i];
    end
    plot_solutions(t, x, h);
end

function y = f(x, k, a, b)
    y = k*(a - x)*(b - x);
end

function plot_solutions(t, x, h)
    plot(t, x);
    hold on;
    plot(t, exact_x(t));
    hold off;
    error = norm(x - exact_x(t), 'inf');
    title(...
        sprintf(...
        'x(t) (numerical and exact solutions), h = %0.2f, error = %0.3e', ...
        h, error))
    legend({'Numerical solution', 'Exact solution'}, ...
    'Location', 'southeast')
end

function x = exact_x(t)
    x = 350*(1 - exp(-0.2*t))./(7 - 5*exp(-0.2*t));
end

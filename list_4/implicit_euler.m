function [t, y] = implicit_euler(y_0, t_0, T, N, f)
    % Implicit Euler for function y of one variable
    h = (T - t_0)/N;
    t = [t_0];
    y = [y_0];
    syms x;
    for k = 1:N
        y_k_old = y(end);
        t_k = t_0 + k*h;
        y_k = solve(f(t_k, x)*h + y_k_old - x == 0, x);
        y = [y, eval(y_k)];
        t = [t, t_k];
    end
end


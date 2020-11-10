function [t, y] = explicit_euler(y_0, t_0, T, N, f)
    % Explicit (common) Euler for function y of one variable
    h = (T - t_0)/N;
    t = [t_0];
    y = [y_0];
    for k = 1:N
        y_k_old = y(end);
        t_k_old = t(end);
        y_k = y_k_old + f(t_k_old, y_k_old)*h;
        y = [y, y_k];
        t_k = t_0 + k*h;
        t = [t, t_k];
    end
end


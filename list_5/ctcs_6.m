function [x, t, u] = ctcs_6(dx, dt, T)
    % Solve exercise 6 (wave equation and CTCS).
    % Examples:
        % [x, t, u] = ctcs_6(0.01, 0.01, 1);
        % [x, t, u] = ctcs_6(0.01, 0.0105, 1);
    [x, t, u] = ctcs(1, 0, 1, @f, @g, @p, @q, dx, dt, T);
    plot_u(x, T, u)
end

function [x, t, u] = ctcs(c, a, b, f, g, p, q, dx, dt, T)
    % Calculate the generic CTCS method.
    x = a:dx:b;
    x = x';
    [N_x, ~] = size(x);
    N_x = N_x - 1;
    u_0 = f(x(2:(end-1)));
    
    sigma = c*dt/dx;
    m = N_x - 1;
    A = A_matrix(sigma, m);
    aux = zeros(m, 1);
    aux(1) = p(0);
    aux(end) = q(0);
    u_1 = 1/2*A*u_0 + g(x(2:(end-1)))*dt + sigma^2/2*aux;
    
    t = 0:dt:T;
    t = t';
    [N_t, ~] = size(t);
    N_t = N_t - 1;
    u = zeros(m, N_t + 1);
    u(1:end, 1) = u_0;
    u(1:end, 2) = u_1;
    for i = 3:(N_t+1)
        i = i - 1;
        aux = zeros(m, 1);
        aux(1) = sigma^2*p(i*dt);
        aux(end) = sigma^2*q(i*dt);
        u(1:end, i+1) = A*u(1:end, i) - u(1:end, i-1) + aux;
    end
    u = [p(t)'; u; q(t)'];
end

function plot_u(x, T, u)
    % Plot the function through space.
    plot(x, u(1:end, end))
    hold on
    plot(x, sin(2*pi*x).*(sin(2*pi*T) + cos(2*pi*T)), 's')
    hold off
    title(sprintf('Numeric solution for wave equation when t = %.2f', T));
    xlabel('x') 
    ylabel('u')
    legend({'Numerical solution', 'Exact solution'}, ...
    'Location', 'northeast')
end

function A = A_matrix(sigma, m)
    % Compute the A matrix of recurrence.
    A = eye(m)*(2-2*sigma^2);
    A(1:(end-1), 2:end) = A(1:(end-1), 2:end) + eye(m-1)*sigma^2;
    A(2:end, 1:(end-1)) = A(2:end, 1:(end-1)) + eye(m-1)*sigma^2;
end

function y = f(x)
    % Compute f function.
    y = sin(2*pi*x);
end

function y = g(x)
    y = 2*pi*sin(2*pi*x);
end

function y = p(x)
    y = x - x;
end

function y = q(x)
    y = x - x;
end


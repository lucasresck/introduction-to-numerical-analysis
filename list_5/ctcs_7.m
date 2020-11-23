function [x, t, u] = ctcs_7(dx, dt, a, b, T)
    % Examples:
        % [x, t, u] = ctcs_7(0.04, 0.02, -10, 10, 40);
        % [x, t, u] = ctcs_7(0.04, 0.040008, -10, 10, 40);
    [x, t, u] = ctcs(1, a, b, @f, @g, @p, @q, dx, dt, T);
    plot_u(x, t, u)
end

function [x, t, u] = ctcs(c, a, b, f, g, p, q, dx, dt, T)
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
        u(1:end, i+1) = A*u(1:end, i) - u(1:end, i-1) + aux;
    end
    u = [p(t)'; u; q(t)'];
end

function plot_u(x, t, u)
    [m, n] = size(u);
    c = floor(m/100);
    d = floor(n/50);
    x = x(1:c:m);
    t = t(1:d:n);
    u = u(1:c:m, 1:d:n);
    [t, x] = meshgrid(t', x');
    s = surf(x, t, u);
%     s.EdgeColor = 'none';
end

function A = A_matrix(sigma, m)
    A = eye(m)*(2-2*sigma^2);
    A(1:(end-1), 2:end) = A(1:(end-1), 2:end) + eye(m-1)*sigma^2;
    A(2:end, 1:(end-1)) = A(2:end, 1:(end-1)) + eye(m-1)*sigma^2;
end

function y = f(x)
    y = exp(-x.^2);
end

function y = g(x)
    y = x - x;
end

function y = p(x)
    y = x - x;
end

function y = q(x)
    y = x - x;
end


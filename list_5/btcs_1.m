function [x, t, u] = btcs_1()
    % Solve exercise 1.
    [x, t, u] = btcs(1/pi^2, 0, 1, @f, @a, @b, 0.04, 0.01, 0.5);
    plot_u(x, 0.5, u)
end

function [x, t, u] = btcs(c, a, b, f, fa, fb, dx, dt, T)
    % Calculate the generic BTCS method.
    x = a:dx:b;
    x = x';
    [N_x, ~] = size(x);
    N_x = N_x - 1;
    u_0 = f(x(2:(end-1)));
    
    v = c*dt/dx^2;
    m = N_x - 1;
    A = A_matrix(v, m);
    
    t = 0:dt:T;
    t = t';
    [N_t, ~] = size(t);
    N_t = N_t - 1;
    u = zeros(m, N_t + 1);
    u(1:end, 1) = u_0;
    for i = 2:(N_t+1)
        i = i - 1;
        d = zeros(m, 1);
        d(1) = v*fa(i*dt);
        d(end) = v*fb(i*dt);
        u(1:end, i+1) = linsolve(A, u(1:end, i) + d);
    end
    u = [fa(t)'; u; fb(t)'];
end

function plot_u(x, T, u)
    % Plot the function through space.
    plot(u(:, end))
    hold on
    plot(cos(pi*(x - 1/2))*exp(-T))
    hold off
    title('Numeric solution for EDP when t = 0.5');
    xlabel('x') 
    ylabel('u')
    legend({'Numerical solution', 'Exact solution'}, ...
    'Location', 'northeast')
end

function A = A_matrix(v, m)
    % Compute the A matrix of recurrence.
    A = eye(m)*(1 + 2*v);
    A(1:(end-1), 2:end) = A(1:(end-1), 2:end) + eye(m-1)*(-v);
    A(2:end, 1:(end-1)) = A(2:end, 1:(end-1)) + eye(m-1)*(-v);
end

function y = f(x)
    % Compute f function.
    y = cos(pi*(x - 1/2));
end

function y = a(x)
    y = x - x;
end

function y = b(x)
    y = x - x;
end


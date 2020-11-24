function crank_nicolson_3(c, d, L, T, dt, dx)
    % Solve exercise 3.
    % Examples:
        % crank_nicolson_3(1, 9.5, 1, 10, 0.01, 0.01);
        % crank_nicolson_3(1, 10, 1, 100, 0.01, 0.01);
    [x, t, u] = crank_nicolson(c, d, 0, L, T, dt, dx, @f, @a, @b, @A_matrix, @B_matrix);
    plot_u(x, t, u)
end

function [x, t, u] = crank_nicolson(c, d, a, b, T, dt, dx, f, fa, fb, A_matrix, B_matrix)
    % Calculate the generic Crank-Nicolson method.
    x = a:dx:b;
    x = x';
    [N_x, ~] = size(x);
    N_x = N_x - 1;
    u_0 = f(x(2:(end-1)));
    
    v = c*dt/dx^2;
    m = N_x - 1;
    A = A_matrix(v, d, dt, m);
    B = B_matrix(v, d, dt, m);
    
    t = 0:dt:T;
    t = t';
    [N_t, ~] = size(t);
    N_t = N_t - 1;
    u = zeros(m, N_t + 1);
    u(1:end, 1) = u_0;
    for i = 2:(N_t+1)
        i = i - 1;
        aux = zeros(m, 1);
        aux(1) = v*(fa(i*dt) + fa((i-1)*dt));
        aux(end) = v*(fb(i*dt) + fb((i-1)*dt));
        u(1:end, i+1) = linsolve(A, B*u(1:end, i) + aux);
    end
    u = [fa(t)'; u; fb(t)'];
end

function plot_u(x, t, u)
    % Plot the function through space and time.
    [m, n] = size(u);
    c = floor(m/100);
    d = floor(n/50);
    x = x(1:c:m);
    t = t(1:d:n);
    u = u(1:c:m, 1:d:n);
    [t, x] = meshgrid(t', x');
    surf(x, t, u);
    title('Numeric solution for wave equation');
    xlabel('x') 
    ylabel('t') 
end

function A = A_matrix(v, d, dt, m)
    % Compute the A matrix of recurrence.
    A = eye(m)*(2 + 2*v - d*dt);
    A(1:(end-1), 2:end) = A(1:(end-1), 2:end) + eye(m-1)*(-v);
    A(2:end, 1:(end-1)) = A(2:end, 1:(end-1)) + eye(m-1)*(-v);
end

function B = B_matrix(v, d, dt, m)
    % Compute the B matrix of recurrence.
    B = eye(m)*(2 - 2*v + d*dt);
    B(1:(end-1), 2:end) = B(1:(end-1), 2:end) + eye(m-1)*v;
    B(2:end, 1:(end-1)) = B(2:end, 1:(end-1)) + eye(m-1)*v;
end

function y = f(x)
    % Compute f function.
    y = sin(pi/2*x).^2;
end

function y = a(x)
    y = x - x;
end

function y = b(x)
    y = x - x;
end


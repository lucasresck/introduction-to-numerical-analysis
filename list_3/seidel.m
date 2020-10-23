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
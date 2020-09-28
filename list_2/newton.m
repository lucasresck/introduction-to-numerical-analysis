function x = newton(f, f1, x_0, M, d, tol)
    err = tol + 1;
    factor = M/(2*d);
    x = x_0;
    while err > tol
        x_old = x;
        x = x - f(x)/f1(x);
        err = factor*(x - x_old)^2;
    end
    x = eval(x);
end

function x = edo(x_0, tol)
    x = [x_0];
    diff = 10^3;
    n = 2;
    while diff > tol
        x = [x, x_0];
        for i=1:(n-1)
           x(n) = x(n)*(1+i/(n*n));
        end
        diff = x(n) - x(n-1);
        n = n+1;
    end
    x = -x + x_0*sqrt(exp(1));
    plot(x);
end

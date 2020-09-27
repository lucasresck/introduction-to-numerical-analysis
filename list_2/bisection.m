function x = bisection(f, a, b, tol)
    err = tol + 1;
    while err > tol
        x = (a + b)/2;
        fx = f(x);
        if fx == 0
            return
        else
            if f(a)*fx < 0
                b = x;
            else
                a = x;
            end
        end
        err = b - a;
    end
end

function x = regula_falsi(f, a, b, d, tol)
    err = tol + 1;
    while err > tol
        fa = f(a);
        x = a - ((b-a)*fa)/(f(b) - f(a));
        fx = f(x);
        if fx == 0
            return
        else
            if fa*fx < 0
                b = x;
            else
                a = x;
            end
        end
        err = abs(f(x))/d;
    end
    x = eval(x);
end

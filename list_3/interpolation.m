function [p, table] = interpolation(x, y)
    [n, ~] = size(x);
    table = zeros(n, n+1);
    table(1:end, 1) = x;
    table(1:end, 2) = y;
    for j=3:(n+1)
        for i=1:(-j+6)
           table(i, j) = (table(i+1, j-1) - table(i, j-1))/(table(j-2+i, 1)-table(i, 1));
        end
    end
    syms p(z);
    p(z) = 0;
    for i=0:(n-1)
        syms term(z);
        term(z) = table(1, i+2);
        for j=1:i
            term(z) = term(z)*(z-x(j));
        end
        p(z) = p(z) + term(z);
    end
end

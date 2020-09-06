function x_1 = edo(x_0, N)
    x_1 = zeros(N, 1);
    x_1(1) = x_0;
    for n=2:N
        x_1(n) = x_0;
        for i=1:(n-1)
           x_1(n) = x_1(n)*(1+i/(n*n));
        end
    end
end

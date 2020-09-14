function [x, n, time] = sor(w, m)
    tic
    b = full(b_vector(m));
    m2 = m^2;
    x = zeros(m2, 1);
    err = 1;
    n = 0;
    mod_vector = zeros(m2, 1);
    for i=1:m2
        mod_vector(i) = mod(i, m);
    end
    while err > 10e-6
        for i=1:m2
            sum = 0;
            if i > m
                sum = sum - x(i-m);
            end
            if i <= m2-m
                sum = sum - x(i+m);
            end
            mod_i = mod_vector(i);
            if mod_i ~= 0
                sum = sum - x(i+1);
            end
            if mod_i ~= 1
                sum = sum - x(i-1);
            end
            x(i) = x(i) + w*((b(i)-sum)/4-x(i));
        end
        err = norm(x-2, 'inf');
        n = n + 1;
    end
    time = toc;
end

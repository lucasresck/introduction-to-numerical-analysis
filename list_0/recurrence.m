function x = recurrence(n)
    x = zeros(n, 1);
    x(1) = 1;
    x(2) = 1/7;
    for i=3:n
        x(i) = 22/7*x(i-1) - 3/7*x(i-2);
    end
end

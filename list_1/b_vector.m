function b = b_vector(m)
    b = sparse(m^2, 1);
    b(1:m) = 1;
    b((m^2-m+1):m^2) = 1;
    for i=1:m
        b(i*m) = b(i*m) + 1;
        b(i*m-m+1) = b(i*m-m+1) + 1;
    end
    b = b*2;
end

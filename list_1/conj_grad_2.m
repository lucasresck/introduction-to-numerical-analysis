function x = conj_grad_2(A, b, x)
    [m, ~] = size(A);
    r = b - A*x;
    d = r;
    for k=0:(m-1)
        if norm(r, 'inf') == 0
            return
        end
        Ad = A*d;
        alpha = d'*r/(d'*Ad);
        x = x + alpha*d;
        r = r - alpha*Ad;
        beta = r'*A*d/(d'*Ad);
        d = r - beta*d;
    end
end

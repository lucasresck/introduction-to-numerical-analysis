function x_t = interpolator(M_i, M_ip, t_i, t_ip, x_i, x_ip, n)
    k_1 = (x_ip-x_i)*n - (M_ip - M_i)/(6*n);
    k_2 = x_i - M_i/(6*n^2) - k_1*t_i;
    syms x_t(t);
    x_t(t) = n*M_i/6*(t_ip-t)^3 + n*M_ip/6*(t-t_i)^3 + k_1*t + k_2;
end
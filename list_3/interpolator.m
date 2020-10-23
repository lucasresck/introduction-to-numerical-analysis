function x_t = interpolator(M_i, M_ip, t_i, t_ip, x_i, x_ip, n)
    c_1 = (x_ip-x_i)*n - (M_ip - M_i)/(6*n);
    c_2 = x_i - M_i/(6*n^2) - c_1*t_i;
    syms x_t(t);
    x_t(t) = n*M_i/6*(t_ip-t)^3 + n*M_ip/6*(t-t_i)^3 + c_1*t + c_2;
end
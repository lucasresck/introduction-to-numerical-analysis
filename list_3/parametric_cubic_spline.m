function parametric_cubic_spline(x, y, k_0, k_1)
    x_t = anchored_cubic_spline(x, k_0, k_1);
    y_t = anchored_cubic_spline(y, k_0, k_1);
    plot(x_t, y_t, '-k');
end
function [x, y] = interactive_spline(n, k, k_0, k_1)
    clf;
    axis([0, 10, 0, 10]);
    hold on;
    [x_click, y_click] = ginput(1);
    x = [x_click];
    y = [y_click];
    plot(x, y, 'o');
    for j = 1:k
        for i = 1:n+2
            [x_click, y_click] = ginput(1);
            x = [x; x_click];
            y = [y; y_click];
            plot(x, y, 'o');
        end
        parametric_cubic_spline(x, y, k_0, k_1);
        x = [x(end-1)];
        y = [y(end-1)];
    end
    hold off;
end
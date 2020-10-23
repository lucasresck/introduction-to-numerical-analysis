function interactive_spline(n, k, k_0, k_1)
    clf;
    axis([0, 10, 0, 10]);
    hold on;
    [x_click, y_click] = ginput(1);
    x = [x_click];
    y = [y_click];
    plot(x, y, 'o', ...
        'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'k');
    for j = 1:k
        for i = 1:n+2
            [x_click, y_click] = ginput(1);
            x = [x; x_click];
            y = [y; y_click];
            if i == 1 || i == n+2
                plot(x_click, y_click, 'o', ...
                'MarkerEdgeColor', 'r', ...
                'MarkerFaceColor', 'r');
                plot(x(end-1:end), y(end-1:end), '--', ...
                    'Color', [200, 200, 200]/255);
            else
                if i == n+1
                    plot(x_click, y_click, 'o', ...
                    'MarkerEdgeColor', 'k', ...
                    'MarkerFaceColor', 'k');
                else
                    plot(x_click, y_click, 'o', ...
                    'MarkerEdgeColor', 'b', ...
                    'MarkerFaceColor', 'b');
                end
            end
        end
        parametric_cubic_spline(x, y, k_0, k_1);
        x = [x(end-1)];
        y = [y(end-1)];
    end
    hold off;
end
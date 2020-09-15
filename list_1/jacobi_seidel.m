function [x_jacobi, x_seidel, k_jacobi, k_seidel] = jacobi_seidel(n, epsilon)
    tic
    if n < 4
        error('n is less than 4.');
    else
        if mod(n, 2) ~= 0
            error('n is not an ever number.');
        end
    end
    x_jacobi = zeros(n, 1);
    x_seidel = zeros(n, 1);    
    err_jacobi = 1;
    err_seidel = 1;
    k_jacobi = 0;
    k_seidel = 0;
    while err_jacobi >= epsilon || err_seidel >= epsilon
        % Jacobi needs the old values when iterating
        x_jacobi_old = x_jacobi;
        for i=1:n
            % Here we only continue to iterate Seidel if
            % it's necessary
            if err_seidel >= epsilon
                sum_s = 0;
                if i == n/2 || i == n/2+1
                    sum_s = sum_s + x_seidel(i-1);
                    sum_s = sum_s + x_seidel(i+1) + 1;
                    sum_s = sum_s/3;
                else
                    if i == 1
                        sum_s = sum_s + x_seidel(2);
                        sum_s = sum_s - 1/2*x_seidel(n) + 5/2;
                        sum_s = sum_s/3;
                    else
                        if i == n
                            sum_s = sum_s + x_seidel(n-1);
                            sum_s = sum_s - 1/2*x_seidel(1) + 5/2;
                            sum_s = sum_s/3;
                        else
                            sum_s = sum_s + x_seidel(i-1);
                            sum_s = sum_s + x_seidel(i+1);
                            sum_s = sum_s - 1/2*x_seidel(n+1-i) + 3/2;
                            sum_s = sum_s/3;
                        end
                    end
                end
                x_seidel(i) = sum_s;
            end
            if err_jacobi >= epsilon
                sum_j= 0;
                if i == n/2 || i == n/2+1
                    sum_j = sum_j + x_jacobi_old(i-1);
                    sum_j = sum_j + x_jacobi_old(i+1) + 1;
                    sum_j = sum_j/3;
                else
                    if i == 1
                        sum_j = sum_j + x_jacobi_old(2);
                        sum_j = sum_j - 1/2*x_jacobi_old(n) + 5/2;
                        sum_j = sum_j/3;
                    else
                        if i == n
                            sum_j = sum_j + x_jacobi_old(n-1);
                            sum_j = sum_j - 1/2*x_jacobi_old(1) + 5/2;
                            sum_j = sum_j/3;
                        else
                            sum_j = sum_j + x_jacobi_old(i-1);
                            sum_j = sum_j + x_jacobi_old(i+1);
                            sum_j = sum_j - 1/2*x_jacobi_old(n+1-i) + 3/2;
                            sum_j = sum_j/3;
                        end
                    end
                end
                x_jacobi(i) = sum_j;
            end
        end
        % We do not update k nor time nor error for Jacobi if
        % it's not necessary anymore
        if err_jacobi >= epsilon
            err_jacobi = norm(x_jacobi-1,'inf');
            k_jacobi = k_jacobi + 1;
            toc_jacobi = toc;
        end
        % The same for Seidel
        if err_seidel >= epsilon
            err_seidel = norm(x_seidel-1,'inf');
            k_seidel = k_seidel + 1;
            toc_seidel = toc;
        end
    end
    display([toc_jacobi, toc_seidel]);
end
function [x_jacobi, x_seidel, k_jacobi, k_seidel] = jacobi_seidel(n, epsilon)
    if n <= 0
        error('n is not positive.');
    else
        if mod(n, 2) ~= 0
            error('n is not an ever number.');
        end
    end
    x_jacobi = zeros(n, 1);
    x_seidel= zeros(n, 1);
    
    err_jacobi = 1;
    err_seidel = 1;
    k_jacobi = 0;
    k_seidel = 0;
    while err_jacobi >= epsilon || err_seidel >= epsilon
        x_jacobi_old = x_jacobi;
        for i=1:n
            if err_seidel >= epsilon
                sum_seidel = 0;
                if i == n/2 || i == n/2+1
                    sum_seidel = sum_seidel + x_seidel(i-1) + x_seidel(i+1) + 1;
                    sum_seidel = sum_seidel/3;
                else
                    if i == 1
                        sum_seidel = sum_seidel + x_seidel(2) - 1/2*x_seidel(n) + 5/2;
                        sum_seidel = sum_seidel/3;
                    else
                        if i == n
                            sum_seidel = sum_seidel + x_seidel(n-1) - 1/2*x_seidel(1) + 5/2;
                            sum_seidel = sum_seidel/3;
                        else
                            sum_seidel = sum_seidel + x_seidel(i-1) + x_seidel(i+1) - 1/2*x_seidel(n+1-i) + 3/2;
                            sum_seidel = sum_seidel/3;
                        end
                    end
                end
                x_seidel(i) = sum_seidel;
            end
            if err_jacobi >= epsilon
                sum_jacobi= 0;
                if i == n/2 || i == n/2+1
                    sum_jacobi = sum_jacobi + x_jacobi_old(i-1) + x_jacobi_old(i+1) + 1;
                    sum_jacobi = sum_jacobi/3;
                else
                    if i == 1
                        sum_jacobi = sum_jacobi + x_jacobi_old(2) - 1/2*x_jacobi_old(n) + 5/2;
                        sum_jacobi = sum_jacobi/3;
                    else
                        if i == n
                            sum_jacobi = sum_jacobi + x_jacobi_old(n-1) - 1/2*x_jacobi_old(1) + 5/2;
                            sum_jacobi = sum_jacobi/3;
                        else
                            sum_jacobi = sum_jacobi + x_jacobi_old(i-1) + x_jacobi_old(i+1) - 1/2*x_jacobi_old(n+1-i) + 3/2;
                            sum_jacobi = sum_jacobi/3;
                        end
                    end
                end
                x_jacobi(i) = sum_jacobi;
            end
        end
        if err_jacobi >= epsilon
            err_jacobi= norm(x_seidel-1,'inf');
            k_jacobi = k_jacobi + 1;
        end
        if err_seidel >= epsilon
            err_seidel = norm(x_seidel-1,'inf');
            k_seidel = k_seidel + 1;
        end
    end
end

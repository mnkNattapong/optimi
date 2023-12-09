[xmin, fmin, Xk, Fk, Gk, Lk, nF, nG, IFLAG] = BFGS(@Rosenbrock, [2;5], 1e-6, 1e-4, 0.98, 1000);

nF = [0 nF];   % there is no cumputation of f at iteration 0th
nG = [0 nG];   % there is no cumputation of f at iteration 0th

fprintf('% 5s % 13s % 13s % 15s % 15s % 15s % 15s % 13s % 13s \n', ...
    'Iter', 'x1', 'x2', 'f', 'grad_f1', 'grad_f2', 'lambda', 'nF', 'nG')

for i = 1:length(Xk)
    fprintf('% 5.2d % 13.7f % 13.7f % 15.5f % 15.5f % 15.5f % 15.5f % 13.5f % 13.5f \n', ...
        i-1, Xk(1,i), Xk(2,i), Fk(i), Gk(1,i), Gk(2,i), Lk(i), nF(i), nG(i))
end

fprintf('total number of f computation is %i \n', sum(nF))
fprintf('total number of grad_f computation is %i \n', sum(nG))


[xmin,fmin,Xk,Fk,Gk,Lk,nF,nG,IFLAG,nReset] = CG(@Rosenbrock, [2;5], 1e-6, 1e-4, 0.25, 1000, 2);

nF = [0 nF];   % there is no cumputation of f at iteration 0th
nG = [0 nG];   % there is no cumputation of f at iteration 0th
nReset = [0 nReset];

fprintf('% 5s % 13s % 13s % 15s % 15s % 15s % 13s % 13s \n', ...
    'Iter', 'x1', 'x2', 'f', 'lambda', 'nF', 'nG', 'nReset')

for i = 1:length(Xk)
    fprintf('% 5.2d % 13.7f % 13.7f % 15.5f % 15.5f % 15.5f % 13.5f % 13.5f \n', ...
        i-1, Xk(1,i), Xk(2,i), Fk(i), Lk(i), nF(i), nG(i), nReset(i))
end

fprintf('total number of f computation is %i \n', sum(nF))
fprintf('total number of grad_f computation is %i \n', sum(nG))

plot(nReset, '-o')
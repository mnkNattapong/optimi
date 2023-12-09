[xmin,fmin,Xk,Fk,Gk,nF,nG,nH,IFLAG] = Newton(@FunctionName,[-1.2;1],1e-6,1e-2,1e-4,100);

% show the result for each iteration
fprintf('% 6s % 10s % 14s % 14s % 21s % 17s \n', 'Iter', 'x_1', 'x_2', 'f', 'gradient_fx1', 'gradient_fx2');
for i = 0:length(Xk)-1
    fprintf('% 5.2d % 14.7f % 14.7f % 16.5f % 16.5f % 16.5f \n', i, Xk(1,i+1), Xk(2,i+1), Fk(i+1), Gk(1,i+1), Gk(2,i+1));
end

% show the final result
disp(['Counts computation of f: ', num2str(nF)])
disp(['Counts computation of gradient of f: ', num2str(nG)])
disp(['Counts computation of Hessian of f: ', num2str(nH)])

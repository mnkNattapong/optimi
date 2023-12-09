% input for the funciton
a = 0;              % lower bound of the interval
b = 0.1;            % upper bound of the interval
epsilon = 1e-10;    % the error tolerance 
itmax = 100;        % maximum iteration tolerance
    
b0 = b;             % record initial input of b to show after the iteration is done

% compute derivative of function at point a
diff_f_a = diff_f(a);

% check the condition if the optimal solution is in the interval [a,b]

% when the optimal solution is not in the interval, 
% the derivative at point b has the same sign as that of at the point a,
% so we increase b until [a,b] contains the optimal solution
while diff_f_a * diff_f(b) > 0      % not contains optimal soltution => same sign => add b
    b = b + 0.01;
end

% when [a,b] contains the optimal solution => call golden
[xmin, fmin, IFLAG, IFunc, Ak, Bk, X1k, X2k] = golden(a, b, epsilon, itmax);

% show the result for each iteration
fprintf('% 6s % 14s % 21s % 20s % 22s \n', 'Iter', 'a', 'x_1', 'x_2', 'b');
for i = 0:length(Ak)-1
    fprintf('% 5.2d % 20.10f % 20.10f % 20.10f % 20.10f \n', i, Ak(i+1), X1k(i+1), X2k(i+1), Bk(i+1));
end

% show the final result
disp(['Input: The interval [a,b] = [', num2str(a), ' ', num2str(b0), ']'] )
disp(['Result: The optimal solution is at x = ', num2str(xmin), ...
       ' which gives the optimal value f(x) = ', num2str(fmin), ...
       ' with IFunc = ', num2str(IFunc)])


% input for the funciton
a = 0;              % lower bound of the interval
b = 0.1;            % upper bound of the interval
epsilon = 1e-10;    % the error tolerance 
itmax = 100;        % maximum iteration tolerance
total_IFunc = 0;    % initial IFunc

b0 = b;             % record initial input of b to show after the iteration is done

% Try the first golden 
% If the optimal solution is in [a,b], so we compute it
[xmin, fmin, IFLAG, IFunc, Ak, Bk, X1k, X2k] = golden(a, b, epsilon, itmax);
total_IFunc = total_IFunc + IFunc;   % IFunc is added from the fist golden

% We will not change the point a, so we calculate f(a) and keep it
f_a = f(a);

% After the first call of golden function, we can compare fmin and f(b)
% If fmin > f(b), it means that xmin between [a,b] is not the optimal solution

% So we assign xmin as the varible x
% Then we increase b just a little bit and see if x which belongs to [a,b] gives
% f(x) that f(x) > f(b) or not

% If f(x) < f(b), we can ensure that the interval [a,b] contains the
% optimal solution

x = xmin;
if fmin > min(f_a, f(b))
    while f(x) > min(f_a, f(b))   % we calculate f(x) and compare to (f_a and f(b))
        b = b + 0.01;             % increase b 
        x = b - epsilon;          % x is the point belongs to [a,b]
    end

    % when [a,b] contains the optimal solution => call golden
    [xmin, fmin, IFLAG, IFunc, Ak, Bk, X1k, X2k] = golden(a, b, epsilon, itmax);
    total_IFunc = total_IFunc + IFunc;    % IFunc is added from the second golden
end

IFunc = total_IFunc    % just change the name of the variable

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

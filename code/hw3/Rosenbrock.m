function [f,gradient] = Rosenbrock(x, options)
% Let x be a point of x of each itearation.
% Define the function of f and g from hand calculation
% for each option we substitute the point x in to each function follows the chosen function
    
    fun = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
    grad_fun = @(x) [-400*x(1)*(x(2)-x(1)^2)-2*(1-x(1));
                                200*(x(2)-x(1)^2)      ];

    if options == 1
        f = fun(x);
        gradient = 0;  % No gradient computation for this option

    elseif options == 2
        f = fun(x);
        gradient = grad_fun(x);

    else
        error('Invalid options. Choose either 1 or 2');

    end

end



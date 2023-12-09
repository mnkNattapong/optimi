function [f, gradient, Hessian] = FunctionName(x, options)
% Let x be a point of x of each itearation.
% Define the function of f, g, and H from hand calculation
% for each option we substitute the point x in to each function follows the
% chosen function
        fun = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
        grad_fun = @(x) [-400*x(1)*(x(2)-x(1)^2)-2*(1-x(1));
                                   200*(x(2)-x(1)^2)    ];
        hessian_fun = @(x) [-400*x(2)+1200*x(1)^2+2  -400*x(1);
                                 -400*x(1)              200   ];
    if options == 1
        f = fun(x);
        gradient = 0;  % No gradient computation for this option
        Hessian = 0;   % No Hessian computation for this option

    elseif options == 2
        f = fun(x);
        gradient = grad_fun(x);
        Hessian = 0;  % No Hessian computation for this option

    elseif options == 3
        f = fun(x);
        gradient = grad_fun(x);
        Hessian = hessian_fun(x);

    else
        error('Invalid options. Choose from 1, 2, or 3.');
        
    end
end

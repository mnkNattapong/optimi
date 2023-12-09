function [xmin,fmin,IFLAG,nF,nG] = golden(FunctionName,a,b,ep_rel,ep_abs,s,itmax)
    T = (sqrt(5)-1)/2;          % golden ratio
    k = 0;                      % initial number of iterations
    
    x1 = a + (1-T)*(b-a);       % initial x1
    x2 = b - (1-T)*(b-a);       % initial x2

    [f_a, grad_fa, hess_fa] = FunctionName(a,1);           
    [f_b, grad_fb, hess_fb] = FunctionName(b,2);            

    [f_x1, grad_fx1, hess_fx1] = FunctionName(x1,2);
    [f_x2, grad_fx2, hess_fx2] = FunctionName(x2,2);

    nF = 4;  nG = 3;            % initial nF and nG from the computation above
    
    if dot(s,grad_fb) > 0        % check if xmin belongs to [a,b]

        % while loop to find the optimal solution until the condition is satisfied
        while (abs(dot(s,grad_fx2)) > abs(dot(s,grad_fx1))*ep_rel + ep_abs) && k < itmax
            k = k + 1;
            
            % similar to previous homework, follows the algorithm
            if f_x2 > f_x1
                b = x2;
            else
                a = x1;
            end
            
            x1 = a + (1-T)*(b-a);
            x2 = b - (1-T)*(b-a);
            
            % compute f at x1 and x2
            [f_x1, grad_fx1, hess_fx1] = FunctionName(x1,2);
            [f_x2, grad_fx2, hess_fx2] = FunctionName(x2,2);
            
            % number of nF and nG are increased from the above computation
            nF = nF + 2;
            nG = nG + 2;
        end   

        % after iteration is done (raise out the FLAG)
        if k == itmax     
            IFLAG = -999;
            disp('golden iterations exceed')
       
        else
            IFLAG = 0;
            disp('golden optimum is found')
        end

        xmin = (x1+x2)/2;                      % average optimal solution by the middle point
        [fmin, grad_fmin, hess_fmin] = FunctionName(xmin,1);   % evaluate optimal value
        nF = nF + 1;                           % number of computation is added by 1

    else % when xmin does not belong to [a,b]
        IFLAG = -999;
        xmin = 0;
        fmin = 0;
    end
end
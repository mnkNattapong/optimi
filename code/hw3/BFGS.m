function  [xmin,fmin,Xk,Fk,Gk,Lk,nF,nG,IFLAG] = BFGS(FcnName,x0,epsilon,mu,eta,itmax)

    Xk = []; % the array containing iterates Xk
    Fk = []; % the array containing f(xk)
    Gk = []; % the array containing grad_f(xk)
    Lk = []; % the array containing lambda(xk)
    nF = []; % the array containing the number of evaluation of f
    nG = []; % the array containing the number of evaluation of grad_f

    IFLAG = -999; % IFLAG: indicate the success
    
    B = eye(2); % Initial B is set to be an identity matrix
    
    for i = 1:itmax
        num_F = 0;  % the number of evaluation of f
        num_G = 0;  % the number of evaluation of grad_f

        [f0, g0] = FcnName(x0, 2); 
        num_F = num_F + 1; 
        num_G = num_G + 1;
    
        a = 1;       % this is the first given value of lambda
        s = B\(-g0); % set line-search direction
        
        % find lambda satisfied both Armijo's and Strong Wolfe's condition
        [lambda,num_Flin,num_Glin] = linesearch(FcnName, x0, s, a, mu, eta);
    
        num_F = num_F + num_Flin; 
        num_G = num_G + num_Glin;
    
        % store values
        Xk(:,i) = x0; 
        Fk(i) = f0; 
        Gk(:,i) = g0; 
        Lk(i) = lambda;
       
        % update B 
        x1 = x0 + lambda*s;
        [f1,g1] = FcnName(x1, 2); 
        num_F = num_F + 1; 
        num_G = num_G + 1;
    
        diff_g = g1 - g0;
        diff_x = lambda*s;
    
        B = B + diff_g*diff_g'/dot(diff_g,diff_x) - B*(diff_x*diff_x')*B/(diff_x'*B*diff_x);
    
        % set the terminal condition
        if norm(g1) < epsilon % at local minimum, gradient converges to 0
            xmin = x1; 
            fmin = f1; 
            IFLAG = 0;
            disp('search successful.');
            break
        end
    
        % update values for the next iteration
        x0 = x1; 
        f0 = f1; 
        g0 = g1;
        
        nF(i) = num_F;
        nG(i) = num_G;
    end
    
    if IFLAG == -999
        xmin = 0;
        fmin = 0;
        disp('search unsuccessful.');
    end

end




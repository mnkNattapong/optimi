function  [xmin,fmin,Xk,Fk,Gk,Lk,nF,nG,IFLAG,nReset] = CG(FcnName,x0,epsilon,mu,eta,itmax,option)
    Xk = []; % the array containing iterates Xk
    Fk = []; % the array containing f(xk)
    Gk = []; % the array containing grad_f(xk)
    Lk = []; % the array containing lambda(xk)
    nF = []; % the array containing the number of evaluation of f
    nG = []; % the array containing the number of evaluation of grad_f
   
    IFLAG = -999; % IFLAG: indicate the success
    
    nReset = []; % the array containing the record of resetting sk
    
    num_F = 0;
    num_G = 0;

    [f0, g0] = FcnName(x0, 2);
    num_F = num_F + 1;
    num_G = num_G + 1;
    
    % set the initial search direction
    s = -g0;

    for i = 1:itmax
    
        % Perform Line Search to obtain lambda
        a = 1;
        % find lambda satisfied both Armijo's and Strong Wolfe's condition
        [lambda,num_Flin,num_Glin] = linesearch(FcnName, x0, s, a, mu, eta);
        num_F = num_F + num_Flin; 
        num_G = num_G + num_Glin;
        
        % update value x
        x1 = x0 + lambda*s;

        [f1, g1] = FcnName(x1, 2);
        num_F = num_F + 1;
        num_G = num_G + 1;
        
            
        if option == 1      % FR
            B = norm(g1)/norm(g0);
        elseif option == 2  % PR
            B = g1'*(g1-g0)/norm(g0);
        else
            disp('invalid option')
            IFLAG = -999;
            break
        end
        
        % update search direction
        s = -g1 + B*s;

        % Reset criteria
        if dot(s,g1) >= 0    % search direction and gradient is in the same direction
            nReset(i) = 2;
            s = -g1;
        elseif acosd((dot(s,-g1)/norm(s)/norm(-g1))) >= 85  % too large angle
            nReset(i) = 1;
            s = -g1;
        else
            nReset(i) = 0;
        end
      
        % store values
        Xk(:,i) = x0; 
        Fk(i) = f0; 
        Gk(:,i) = g0; 
        Lk(i) = lambda;
        nF(i) = num_F; 
        nG(i) = num_G;
    
        % terminate
        if norm(x1-x0) < epsilon % at local minimum, gradient converges to 0.
            xmin = x1; 
            fmin = f1; 
            IFLAG = 0;
            disp('search successful');
            break
        end

        % update values for the next iteration
        x0 = x1; 
        f0 = f1; 
        g0 = g1; 

        num_F = 0;
        num_G = 0;
    end

    if IFLAG == -999
        xmin = 0; 
        fmin = 0; 
        disp('search unsuccessful');
    end
end










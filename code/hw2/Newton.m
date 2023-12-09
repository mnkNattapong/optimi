function [xmin,fmin,Xk,Fk,Gk,nF,nG,nH,IFLAG] = Newton(FunctionName,x0,epsilon,ep_rel,ep_abs,itmax)
    Xk = []; Fk = []; Gk = [];      % arrays containing f values            
    nF = 0; nG = 0; nH = 0;         % initial counts number of calculation
    IFLAG = 0;                      % initial FLAG

    for i = 1:itmax   % start for loop until it exceeds the maximum number of iteration
        [f0, grad_f0, Hessian_f0] = FunctionName(x0,3);  % compute initial f values

        % number of counts is increased from the computation
        nF = nF+1;
        nG = nG+1;
        nH = nH+1;
        
        % compute the search direction
        s = -Hessian_f0\grad_f0;
        
        % see if s has descent property
        if dot(s,grad_f0) > 0  % s does not have descent property
            s = -grad_f0;      % use steepest descent instead of Newton
        end
        
        % next we find lambdha by using Golden Section Method
        % we use for loop to assign value lamdha as 1.5^j (to be possible to cover the minimum)
        % and then find the optimal lamdha (assign as 1.5^j because it has proper length of
        % increasing and decreasing interval)
        for j = -50:50 
            a = x0;
            b = x0 + (1.5^j)*s;
            
            % try to call golden (in the golden, it already 
            % checks in the function whether xmin belongs to [a,b] or not)
            [x1,f1,IFLAG_lin,nF_lin,nG_lin] = golden(@FunctionName,a,b,ep_rel,ep_abs,s,itmax);
            
            % number of counts is increased from the line search
            nF = nF + nF_lin;
            nG = nG + nG_lin;
            if IFLAG_lin ~= -999
                disp('The optimal solution is found.')
                disp(['use j = ' num2str(j)])
                break
            end
        end


        if IFLAG_lin == -999
            disp("It is quite difficult to find the optimal solution"); 
            IFLAG = -999;
        end
        
        % record the f values of each iteration
        Xk(:,i) = x0;
        Fk(:,i) = f0;
        Gk(:,i) = grad_f0;

        % minimum search for the last iteration
        if norm(x1-x0) < epsilon 
            % calculate the last f values
            [f1, grad_f1, Hessian_f1] = FunctionName(x1,2);
            xmin = x1;
            fmin = f1;
            % record the f values of each iteration
            Xk(:, i+1) = x1;
            Fk(:, i+1) = f1;
            Gk(:, i+1) = grad_f1;
            % number of counts is increased
            nF = nF + 1;
            nG = nG + 1;
            break
        end
        x0 = x1;   % just change the varible to be x0
    end
   
    % when Xk is diverged, raised IFLAG -999
    if i == itmax;
        disp('The function is diverged')
        IFLAG = -999;
    end
       
end
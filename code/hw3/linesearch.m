function [lambda,num_Flin,num_Glin] = linesearch(FcnName, x0, s, a, mu, eta)
    num_Flin = 0;   
    num_Glin = 0;
    a_left = 0;    % lambda = 0
    a_right = a;   % lambda = 1 (this lambda is going to be decreased by back tracking)
    a_pre = a;     % lambda = 1 (this lambda is going to be increased somewhat explaining below) 

    % Firstly, we find lambda which is satisfied Armijo's condition
    [f0, g0] = FcnName(x0, 2);   % compute RHS of the inequality

    x1 = x0 + a_right*s; 
    [f1, ~] = FcnName(x1, 1);    % compute LHS of the inequality

    num_Flin = num_Flin + 2;
    num_Glin = num_Glin + 1;
    
    slope_x0 = dot(s,g0);
    % we initiate lambda = 1, if f1 is higher than f0 (not satisfied
    % Armijo's condition), then we use back tracking algorithm
    
    while f1 > f0 + mu*a_right*slope_x0
        a_right = a_right/3;    % back tracking algorithm
        
        x1 = x0 + a_right*s;    % update x and check again if it satisfies Armijo's condition
        [f1, ~] = FcnName(x1, 1);
        num_Flin = num_Flin + 1;
    end
    % right now we obtained lambda satisfied Armijo's condition
    
    % Check if after back tracking lambda satisfies Strong Wolfe's condition or not, 
    % if so => return lambda 

    [f1, g1] = FcnName(x1, 2);
    num_Flin = num_Flin + 1;
    num_Glin = num_Glin + 1;
 
    slope_x1 = dot(s,g1);               % directional derivative at lambda after back tracking
    RHS_StrongWolfe = -eta*slope_x0;    % compute right hand side of Stong Wolfe's condition
    
    % when lambda satisfied Strong Wolfe's condition => return lambda
    if abs(slope_x1) <= RHS_StrongWolfe 
        lambda = a_right;
        return

    % when lambda is not satisfied Strong Wolfe's condition => check the
    % slope (directional derivative) at that lambda after back tracking
    else
        tau = (sqrt(5)-1)/2;
        % when slope is positive, meaning that there exist optimal lambda
        % which also means that there exists lambda which satisfied Strong Wolfe's condition
        % so we implement golden section method to compute lambda
        if slope_x1 > 0
            % Golden between a_left & a_right
            a = a_left;
            b = a_right;
            c = b - tau*(b - a);
            d = a + tau*(b - a);
            [fc, ~] = FcnName(x0+c*s, 1);
            [fd, ~] = FcnName(x0+d*s, 1);
            num_Flin = num_Flin + 2;
   
            while abs(dot(s,g1)) > RHS_StrongWolfe
                if fd > fc
                    a_right = c;
                    x1 = x0+a_right*s;
                    b = d;
                    d = c;
                    c = a + (1-tau)*(b-a);
                else
                    a_right = d;
                    x1 = x0+a_right*s;
                    a = c;
                    c = d;
                    d = b - (1-tau)*(b-a);
                end
                [f1, g1] = FcnName(x1,2);
                num_Flin = num_Flin + 1;
                num_Glin = num_Glin + 1;

                [fc, ~] = FcnName(x0+c*s,1);
                [fd, ~] = FcnName(x0+d*s,1);
                num_Flin = num_Flin + 2;
            end
            lambda = a_right;  % finally we obtained lambda which satisfied both Armijo and Strong Wolfe 
        
        % when slope is negative, meaning that back tracking algorithm
        % decreased too much value of lambda,
        % how ever we can also implement golden section method to compute lambda
        else
            % in case of the slope at a_pre(lambda = 1) is also negative
            % (we extend it by multiple of 2 in order to ensure that 
            % it covers optimal lambda when we implement golden section method 
            [f_pre,g_pre] = FcnName(x0+a_pre*s,2);
            num_Flin = num_Flin + 1;
            num_Glin = num_Glin + 1;
            while dot(s,g_pre) <= 0
                a_pre = 2*a_pre;
                [f_pre,g_pre] = FcnName(x0+a_pre*s,2);
                num_Flin = num_Flin + 1;
                num_Glin = num_Glin + 1;
            end

            % Golden between a_right & a_pre
            a = a_right;
            b = a_pre;
            c = b - tau*(b - a);
            d = a + tau*(b - a);
            [fc, ~] = FcnName(x0+c*s,  1);
            [fd, ~] = FcnName(x0+d*s, 1);
            num_Flin = num_Flin + 2;
            while abs(dot(s,g1)) > RHS_StrongWolfe

                if fd > fc
                    a_pre = c;
                    x1 = x0+a_pre*s;
                    b = d;
                    d = c;
                    c = a + (1-tau)*(b - a); 
                else
                    a_pre = d;
                    x1 = x0+a_pre*s;
                    a = c;
                    c = d;
                    d = b - (1-tau)*(b - a);
                end
                [f1,g1] = FcnName(x1,2);
                num_Flin = num_Flin + 1;
                num_Glin = num_Glin + 1;

                [fc,~] = FcnName(x0+c*s,1);
                [fd,~] = FcnName(x0+d*s,1);
                num_Flin = num_Flin + 2;
                
            end
            lambda = a_pre; % finally we obtained lambda which satisfied both Armijo and Strong Wolfe 
        end
    end
end
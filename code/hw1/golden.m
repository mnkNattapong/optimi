function [xmin, fmin, IFLAG, IFunc, Ak, Bk, X1k, X2k] = golden(a, b, epsilon, itmax)

    T = (sqrt(5)-1)/2;          % golden ratio
    k = 0;                      % initial number of iterations
    
    x1 = a + (1-T)*(b-a);       % initial x1
    x2 = b - (1-T)*(b-a);       % initial x2

    IFunc = 0;                  % initial IFunc
    
    Ak = [];  Ak(1) = a;        % assign all initial data into the vector
    Bk = [];  Bk(1) = b;
    X1k = []; X1k(1) = x1;
    X2k = []; X2k(1) = x2;
    
    % start the iteration 
    while ((abs(b-a) > epsilon) && (k < itmax))
            k = k + 1;
            
            % follows the pseudo code
            if f(x2) > f(x1)       % use a file name 'f.m' to compute f(x)
                 b    =  x2;       
                 x2   =  x1;               % next x2
                 x1   =  a + (1-T)*(b-a);  % next x1
                f_x1  =  f(x1);
                f_x2  =  f(x2);

            else
                 a    =  x1;
                 x1   =  x2;               % next x1
                 x2   =  b - (1-T)*(b-a);  % next x2
                f_x1  =  f(x1);
                f_x2  =  f(x2);

            end
            
            % record the result of each iteration into the vector
            Ak(k+1)  = a;
            Bk(k+1)  = b;
            X1k(k+1) = x1;
            X2k(k+1) = x2;
              IFunc  = IFunc + 2;  % each iteration computes 2 evaluations (f_x1 & f_x2)
    end
    
    xmin = (x1+x2)/2;       % average optimal solution by the middle point
    fmin = f(xmin);         % evaluate optimal value
    IFunc = IFunc + 1;      % number of computation is added by 1 (the last computation)

    % after iteration is done (raise out the FLAG)
    if k == itmax     
        IFLAG = -999;
        disp('number of iterations exceed')
   
    else
        IFLAG = 0;
        disp('the optimum is found')
    
    end

end
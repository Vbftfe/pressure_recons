function pressureField = SORPoissonSolver(param,pivData,sourceTerm,currentTime)


    %% SETTINGS
    ux = pivData.ux{currentTime};
    uy = pivData.ux{currentTime};
    
    % Length of array representing 2D grid
    arrayLength = length(pivData.x);
    
    % x and y width of 2D grid
    Lx = length(unique(pivData.x));
    Ly = arrayLength/Lx;
    
    % Non-dimensional steps
    dx = 1; % = dy
    
    % SOR parameters
    maxIter = param.maxIter;                    % Maximum iterations
    minError = param.minError;                  % Convergence criteria
    if strcmp(param.omega,'opt')
        xi = (0.5*(cos(pi/Lx) + cos(pi/Ly)))^2;
        omegaOpt = 2*(1-sqrt(1-xi))/xi;             % Optimal SOR parameter for fastest convergence
        omega = omegaOpt;                           % Choice of SOR parameter
    else
        omega = param.omega;
    end
    
    % Initial pressure field
%     pressureField = zeros(1,arrayLength);
%     pressureField = ones(1,arrayLength)*1.22e04;

    %% BOUNDARY CONDITIONS
    
    % Set Dirichlet boundary conditions (assuming farfield Bernoulli
    % pressure)
    % Obtain domain boundaries
    [domainBoundary,internalDomain] = getDomainBoundary(pivData);
    
%     % Obtain velocity at boundaries
%     uxBC = pivData.ux{currentTime}(domainBoundary);
%     uyBC = pivData.uy{currentTime}(domainBoundary);
%     
%     % Calculate Bernoulli pressure at boundaries
%     pressureBC = 0.5*(1 - (uxBC.^2 + uyBC.^2));
%     pressureField(domainBoundary) = pressureBC;
    % 获取边界处的压力�?
%     pressureField(domainBoundary) = param.turePressure(domainBoundary);
%     pressureField((end - Lx + 1):end) = param.turePressure(1:Lx);
    
    % Convert to 2D grid
    pressureField = array2grid(pivData,pressureField);
%     oldField = pressureField;
%     domainBoundary = array2grid(pivData,domainBoundary);
    internalDomain = array2grid(pivData,internalDomain);
    
%     domain = padarray(domainBoundary + internalDomain, [1, 1]);
%     pressureField = padarray(pressureField, [1, 1]);
    
    % Iterate algorithm towards convergence
    for iter = 1:maxIter
        % Record previous field before calculation
        oldField = pressureField;
    
        % Loop through grid
        for j = 1:Ly
            for i = 1:Lx

                % Bottom left starting point
                % Perform computations inside bounded domain
                if internalDomain(j,i)
%                 if domain(j,i)

                    % Calculate pressure change
                    dp = 0.25*(pressureField(j,i-1) + pressureField(j,i+1)...
                        + pressureField(j-1,i) + pressureField(j+1,i) - sourceTerm(j,i)*dx^2);
%                     dp = 0.25*(pressureField(j,i-1) + oldField(j,i+1)...
%                         + pressureField(j-1,i) + oldField(j+1,i) - sourceTerm(j,i)*dx^2);

                    % Calculate new pressure field
                    pressureField(j,i) = (1-omega)*pressureField(j,i) + omega*dp;

                end
            end

        end
        
        % Check convergence (take minimum of difference and relative error)
%         error = max(abs(pressureField - oldField)./abs(oldField),[],'all');
        diff_error = max(abs(pressureField - oldField),[],'all');
%         if error < minError %|| diff_error < minError
%             break;
%         end
        if diff_error < minError
            break;
        end
        
         % 将新的边界条件带�?
        boundary_pressure = calc_boundary_pressure2(ux, uy, grid2array(pressureField), domainBoundary, param, Lx, Ly);
        pressureField = grid2array(pressureField);
        pressureField(domainBoundary) = boundary_pressure(domainBoundary);
        pressureField = array2grid(pivData, pressureField);
%         pressureField((end - Lx + 1):end) = param.turePressure(1:Lx);
    end
        
    if iter == maxIter
        disp('Solution has not converged');
    end
    fprintf("totally iterate %d times\n", iter);
    
%     % Return pressure coefficient
%     pressureField = 2*pressureField;
%     pressureField = imcrop(pressureField, [2, 2, Lx - 1, Ly - 1]);

end
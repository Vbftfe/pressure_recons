function [domainBoundary,internalDomain] = myGetDomainBoundary(pivData)
    % Convert domain into 2D grid
    fluidDomain = array2grid(pivData,pivData.domain);
    
    
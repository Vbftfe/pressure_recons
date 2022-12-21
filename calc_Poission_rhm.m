function pressureField = calc_Poission_rhm(pivData, param, currentTime)
% 参考文献: https://sci-hub.se/https://www.sciencedirect.com/science/article/pii/0021999187900088
% omega = duy_dx - dux_dy;
[domainBoundary, internalDomain] = getDomainBoundary(pivData);
ux = pivData.ux{currentTime};
uy = pivData.ux{currentTime};
Lx = length(unique(pivData.x));
omega = zeros(length(ux), 1);
velSource = omega;
h = 1;
dx = h;
dy = dx;
% 计算边界上的导数值
[a, dux_dy, duy_dx, d] = calc_boundary_diff(pivData, Lx, domainBoundary, currentTime);
omega(domainBoundary) = duy_dx(domainBoundary) - dux_dy(domainBoundary);

for i = 1:length(ux)
    if internalDomain(i)
        duy_dx = (uy(i + 1) - uy(i - 1))/(2*dx);
        dux_dy = (ux(i + Lx) - uy(i - Lx))/(2*dy);
        omega(i) = duy_dx - dux_dy;
    end
end

for i = 1:length(ux)
    if internalDomain(i)
        velSource(i) = ( (uy(i + 1) + uy(i))*(omega(i + 1) + omega(i)) ...
                       - (uy(i) + uy(i - 1))*(omega(i) + omega(i - 1))...
                       - (ux(i + Lx) + ux(i))*(omega(i + Lx) + omega(i))...
                       + (ux(i) + ux(i - Lx))*(omega(i) + omega(i - Lx)) )*h/4;
    end
end

pressureField = SORPoissonSolver(param,pivData,velSource,currentTime);
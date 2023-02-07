function [pressure_gd_x, pressure_gd_y] = calc_pressure_gd(pivData, Lx, Ly, param)
[domainBoundary, internalDomain] = getDomainBoundary(pivData);
ux = pivData.ux{currentTime};
uy = pivData.uy{currentTime};
miu = param.miu;
rho = param.rho;
dx = 1;
dy = dx;
dux_dx = zeros(1,length(ux));
dux_dy = dux_dx;
duy_dx = dux_dx;
duy_dy = dux_dx;
dux_dx_2 = dux_dx;
dux_dy_2 = dux_dx;
duy_dx_2 = dux_dx;
duy_dy_2 = dux_dx;
pressure_gd_x = dux_dx;
pressure_gd_y = dux_dx;

for i = 1:length(ux)
    if internalDomain(i)
        % Compute first derivatives (central difference scheme)
        dux_dx(i) = (ux(i+1)-ux(i-1))/(2*dx);
        dux_dy(i) = (ux(i+Lx)-ux(i-Lx))/(2*dy);
        duy_dx(i) = (uy(i+1)-uy(i-1))/(2*dx);
        duy_dy(i) = (uy(i+Lx)-uy(i-Lx))/(2*dy);
        dux_dx_2(i) = ux(i + 2) + ux(i) - 2*ux(i + 1);
        duy_dx_2(i) = uy(i + 2) + uy(i) - 2*uy(i + 1);
        dux_dy_2(i) = ux(i + 2*Lx) + ux(i) - 2*ux(i + Lx);
        duy_dy_2(i) = uy(i + 2*Lx) + uy(i) - 2*uy(i + Lx);
    end
end

pressure_gd_x(internalDomain) = miu*(dux_dx_2 + dux_dy_2)...
    - rho*(ux.*dux_dx + uy.*dux_dy);
pressure_gd_y(internalDomain) = miu*(duy_dx_2 + duy_dy_2)...
    - rho*(ux.*duy_dx + uy.*duy_dy);
% 处理边界
[boundary_pressure_gd_x, boundary_pressure_gd_y] = calc_boundary_pressure_gd(ux, uy, domain_boundary, miu, Lx, Ly);
pressure_gd_x(domainBoundary) = boundary_pressure_gd_x;
pressure_gd_y(domainBoundary) = boundary_pressure_gd_y;
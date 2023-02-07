function boundary_pressure = calc_boundary_pressure2(ux, uy, pressure_field, domain_boundary, param, Lx, Ly)
omega = param.vorticity;
% miu = param.miu;
Re = param.Re0;
boundary_pressure = zeros(1, length(ux));
h = 1;

for i = 1:length(ux)
    if domain_boundary(i)
        %% 边界上的计算导数
        % 计算上下边界之外的左边界
        if mod(i - 1, Lx) < Lx/2 ...
            && ~(floor((i - 1) / Lx) == 0 ...
            || floor((i - 1) / Lx) == Ly - 1)
            
            diff = -(omega(i + 1 + Lx) + omega(i + Lx) - omega(i + 1 - Lx) ...
                - omega(i - Lx))/(4*Re) + uy(i + 1)*(omega(i + 1) + omega(i))*h/4;
            boundary_pressure(i) = pressure_field(i + 1) - diff;
        end
        
        % 计算上下边界之外的右边界
        if mod(i - 1, Lx) > Lx/2 ...
            && ~(floor((i - 1) / Lx) == 0 ...
            || floor((i - 1) / Lx) == Ly - 1)
            
            diff = (omega(i + Lx) + omega(i - 1 + Lx) - omega(i - Lx) ...
                - omega(i - 1 - Lx))/(4*Re) - uy(i - 1)*(omega(i) + omega(i - 1))*h/4;
            boundary_pressure(i) = pressure_field(i - 1) - diff;
        end
        
        % 计算左右边界之外的下边界
        if floor((i - 1)/Lx) == 0 && (mod(i, Lx) ~= 0 && mod(i, Lx) ~= 1)
            diff = (omega(i + 1 + Lx) + omega(i + 1) - omega(i - 1 + Lx) ...
                - omega(i - 1))/(4*Re) - ux(i + Lx)*(omega(i + Lx) + omega(i))*h/4;
            boundary_pressure(i) = pressure_field(i + Lx) - diff;
        end
        
        % 计算左右边界之外的上边界
        if floor((i - 1)/Lx) == Ly - 1 && (mod(i, Lx) ~= 0 && mod(i, Lx) ~= 1)
            diff = -(omega(i + 1) + omega(i + 1 - Lx) - omega(i - 1) ...
                - omega(i - 1 - Lx))/(4*Re) + (ux(i) + ux(i - Lx))*(omega(i) + omega(i - Lx))*h/4;
            boundary_pressure(i) = pressure_field(i - Lx) - diff;
        end
    end
end
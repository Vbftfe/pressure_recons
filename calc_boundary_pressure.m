function boundary_pressure = calc_boundary_pressure(ux, uy, pressure_field, domain_boundary, miu, Lx, Ly)
dux_dx_2 = zeros(1, length(ux));
dux_dy_2 = dux_dx_2;
duy_dx_2 = dux_dx_2;
duy_dy_2 = dux_dx_2;
boundary_pressure = dux_dx_2;

for i = 1:length(ux)
    if domain_boundary(i)
%         x = mod(i - 1, Lx);
%         y = floor(i - 1, Lx);
        %% 边界上的计算导数
        % 计算左右边界上x方向的导数
        if mod(i - 1, Lx) < Lx / 2  % 处理左边界
            dux_dx_2(i) = ux(i + 2) + ux(i) - 2*ux(i + 1);
            duy_dx_2(i) = uy(i + 2) + uy(i) - 2*uy(i + 1);
        else % 处理右边界
            dux_dx_2(i) = ux(i) + ux(i - 2) - 2*ux(i - 1);
            duy_dx_2(i) = uy(i) + uy(i - 2) - 2*uy(i - 1);
        end
        
        % 对于上下边界特殊处理
        if floor((i - 1) / Lx) == 0 || floor((i - 1) / Lx) == Ly - 1
            if floor((i - 1) / Lx) == 0 
                % 处理上边界x方向导数(除左上角和右上角)
                if (mod(i, Lx) ~= 1 &&  mod(i, Lx) ~= 0)
                    dux_dx_2(i) = ux(i + 1) + ux(i - 1) - 2*ux(i);
                    duy_dx_2(i) = uy(i + 1) + uy(i - 1) - 2*uy(i);
                end
                % 计算上边界y方向的导数
                dux_dy_2(i) = ux(i) + ux(i + 2*Lx) - 2*ux(i + Lx);
                duy_dy_2(i) = uy(i) + uy(i + 2*Lx) - 2*uy(i + Lx);
            end
            
            if floor((i - 1) / Lx) == Ly - 1 
                % 处理下边界x方向导数(除左下角和右下角)
                if (mod(i, Lx) ~= 1 &&  mod(i, Lx) ~= 0)
                    dux_dx_2(i) = ux(i + 1) + ux(i - 1) - 2*ux(i);
                    duy_dx_2(i) = uy(i + 1) + uy(i - 1) - 2*uy(i);
                end
                % 计算下边界y方向的导数
                dux_dy_2(i) = ux(i) + ux(i - 2*Lx) - 2*ux(i - Lx);
                duy_dy_2(i) = uy(i) + uy(i - 2*Lx) - 2*uy(i - Lx);
            end
            
        else % 除上下边界外，计算边界上关于y方向的导数
            % 对边界为曲面的位置进行特殊处理
            if mod(i - 1, Lx) < Lx / 2
                if isnan(ux(i - Lx))
                    idx = find(ux((i - Lx):end), 1);
                    dux_dy_2(i) = ux(i + Lx) + ux(idx) - 2*ux(i);
                    duy_dy_2(i) = uy(i + Lx) + uy(idx) - 2*uy(i);
                elseif isnan(ux(i + Lx))
                    idx = find(ux((i + Lx):end), 1);
                    dux_dy_2(i) = ux(idx) + ux(i - Lx) - 2*ux(i);
                    duy_dy_2(i) = uy(idx) + uy(i - Lx) - 2*uy(i);
                else
                    dux_dy_2(i) = ux(i + Lx) + ux(i - Lx) - 2*ux(i);
                    duy_dy_2(i) = uy(i + Lx) + uy(i - Lx) - 2*uy(i);
                end
            else
                if isnan(ux(i - Lx))
                    idx = find(ux(1:(i - Lx)), 1);
                    dux_dy_2(i) = ux(i + Lx) + ux(idx) - 2*ux(i);
                    duy_dy_2(i) = uy(i + Lx) + uy(idx) - 2*uy(i);
                elseif isnan(ux(i + Lx))
                    idx = find(ux(1:(i + Lx)), 1);
                    dux_dy_2(i) = ux(idx) + ux(i - Lx) - 2*ux(i);
                    duy_dy_2(i) = uy(idx) + uy(i - Lx) - 2*uy(i);
                else
                    dux_dy_2(i) = ux(i + Lx) + ux(i - Lx) - 2*ux(i);
                    duy_dy_2(i) = uy(i + Lx) + uy(i - Lx) - 2*uy(i);
                end
            end
        end
        
        %% 计算边界的压力值
        % 使用x方向的梯度分量计算左右边界
        if mod(i - 1, Lx) < Lx / 2   % 计算左边界
            boundary_pressure(i) = pressure_field(i + 1) - miu*(dux_dx_2(i) + dux_dy_2(i));
        end
        
        if mod(i - 1, Lx) > Lx / 2   % 计算右边界
            boundary_pressure(i) = pressure_field(i - 1) + miu*(dux_dx_2(i) + dux_dy_2(i));
        end
        
        % 使用y方向的梯度分量计算下边界
        if floor((i - 1) / Lx) == Ly - 1 && (mod(i, Lx) ~= 1 &&  mod(i, Lx) ~= 0)
            boundary_pressure(i) = pressure_field(i - Lx) + miu*(duy_dx_2(i) + duy_dy_2(i));
        end
        
        % 使用y方向的梯度分量计算上边界
        if floor((i - 1) / Lx) == 0 && (mod(i, Lx) ~= 1 &&  mod(i, Lx) ~= 0)
            boundary_pressure(i) = pressure_field(i + Lx) - miu*(duy_dx_2(i) + duy_dy_2(i));
        end
    end
end
% % 将左上角和右上角置零
% boundary_pressure(1) = 0;
% boundary_pressure(Lx) = 0;
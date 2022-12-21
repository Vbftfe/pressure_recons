function [dux_dx, dux_dy, duy_dx, duy_dy] = calc_boundary_diff(pivData, Lx, domainBoundary, currentTime)
    ux = pivData.ux{currentTime};
    uy = pivData.ux{currentTime};
    
    ux(isnan(ux)) = 0;
    uy(isnan(uy)) = 0;
    
    dux_dx = zeros(1,length(ux));
    dux_dy = dux_dx;
    duy_dx = dux_dx;
    duy_dy = dux_dx;
    
    for k = 1:length(ux)
        if domainBoundary(k)
            [dux_dx(k), duy_dx(k)] = calc_du_x(ux, uy, k, Lx);
            [dux_dy(k), duy_dy(k)] = calc_du_y(ux, uy, k, Lx);
        end
    end
end

function [dux_dx, duy_dx] = calc_du_x(ux, uy, idx, Lx)
    if mod(idx - 1, Lx) > Lx / 2
       dux_dx = ux(idx) - ux(idx - 1);
       duy_dx = uy(idx) - uy(idx - 1);
    else
        dux_dx = ux(idx + 1) - ux(idx);
        duy_dx = uy(idx + 1) - uy(idx);
    end
end

function [dux_dy, duy_dy] = calc_du_y(ux, uy, idx, Lx)
    if idx >= 1 && idx <= Lx
       dux_dy = ux(idx + Lx) - ux(idx);
       duy_dy = uy(idx + Lx) - uy(idx);
    else
        dux_dy = ux(idx) - ux(idx - Lx);
        duy_dy = uy(idx) - uy(idx - Lx);
    end
end
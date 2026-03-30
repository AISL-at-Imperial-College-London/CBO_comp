function [gp_J, gp_g1, gp_g2, gp_g3] = fit_all_gps_series(D, lb, ub)
% FIT_ALL_GPS_SERIES Trains surrogate models for the 2-compressor series network.
    
    X = D.U;
    if isempty(X)
        gp_J = []; gp_g1 = []; gp_g2 = []; gp_g3 = []; return;
    end

    lb = lb(:); ub = ub(:); span = ub - lb;
    X_scaled = (X - lb.')./span.'; 

    % Fit 4 individual models
    gp_J  = make_gp(X_scaled, D.cost(:), lb, span, 'se');
    gp_g1 = make_gp_p4(X_scaled, X, D.constr(:,1), lb, span); 
    gp_g2 = make_gp_dts(X_scaled, D.constr(:,2), lb, span, 'se');
    gp_g3 = make_gp_dts(X_scaled, D.constr(:,3), lb, span, 'se');
end
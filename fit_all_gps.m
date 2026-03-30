function [gp_J, gp_g1, gp_g2, gp_g3, gp_g4] = fit_all_gps(D, lb, ub)
% FIT_ALL_GPS Trains all surrogate models for the compressor network.
    
    X = D.U;
    if isempty(X)
        gp_J = []; gp_g1 = []; gp_g2 = []; gp_g3 = []; gp_g4 = [];
        return;
    end

    lb = lb(:);
    ub = ub(:);
    span = ub - lb;
    X_scaled = (X - lb.')./span.'; % Inline scaling

    J_arr  = D.cost(:);
    g1_arr = D.constr(:,1);
    g2_arr = D.constr(:,2);
    g3_arr = D.constr(:,3);
    g4_arr = D.constr(:,4);

    % Fit individual models
    gp_J  = make_gp(X_scaled, J_arr, lb, span, 'se');
    gp_g1 = make_gp_p4(X_scaled, X, g1_arr, lb, span); % X passed as X_orig
    gp_g2 = make_gp_dts(X_scaled, g2_arr, lb, span, 'se');
    gp_g3 = make_gp_dts(X_scaled, g3_arr, lb, span, 'se');
    gp_g4 = make_gp_dts(X_scaled, g4_arr, lb, span, 'se');
end
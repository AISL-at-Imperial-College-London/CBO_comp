function J = objectiveFun_CEI_series(T, P_dis_ref, GP_Power, GP_Pdis, GP_ds1, GP_ds2, Dataset_pressure, Dataset_cost, params)
% OBJECTIVEFUN_CEI_SERIES Computes CEI for the 2-compressor series network.

    T_row = T(:).';   

    % ================= GP Predictions =================
    [mu_J,  std_J]  = gp_predict(GP_Power, T_row);
    [mu_g1, std_g1] = gp_predict(GP_Pdis,  T_row, true); 
    [mu_g2, std_g2] = gp_predict(GP_ds1,   T_row);
    [mu_g3, std_g3] = gp_predict(GP_ds2,   T_row);

    % ================= Expected Improvement (EI) =================
    safe_idxs = Dataset_pressure >= P_dis_ref;
    if ~any(safe_idxs)
        bestJ = 1e7; % Fallback if no feasible points
    else
        bestJ = min(Dataset_cost(safe_idxs));
    end

    sigma = max(std_J, 1e-12);
    improvement = bestJ - mu_J;
    Z = improvement ./ sigma;
    
    EI = improvement .* normcdf_local(Z) + sigma .* normpdf_local(Z);
    EI(sigma < 1e-12) = 0.0;
    EI = max(EI, 0.0);

    % ================= Probability of Feasibility (PoF) =================
    mu_dts  = [mu_g2, mu_g3];
    std_dts = max([std_g2, std_g3], 1e-12);
    
    LCB_dts = mu_dts - params.betaSafe .* std_dts;
    
    PoF_arr = normcdf_local( (LCB_dts - params.surge_safe_margin) ./ std_dts );
    PoF_arr = max(min(PoF_arr, 1.0), 0.0); 
    PoF = prod(PoF_arr); 

    % ================= Tracking Penalty =================
    e = (mu_g1 - P_dis_ref) / params.cei_track_tol;
    P_track_g1 = 1 / (e.^2 + 1);

    % ================= Final Acquisition =================
    J = - (P_track_g1 .* EI .* PoF);
end

% --- Local Math Helpers ---
function y = normpdf_local(x)
    y = 1./sqrt(2*pi) .* exp(-0.5*x.^2);
end
function y = normcdf_local(x)
    y = 0.5 * (1 + erf(x ./ sqrt(2)));
end
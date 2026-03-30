function J = objectiveFun_CEI(T, massflow_ref, GP_Power, GP_Pdis, GP_ds1, GP_ds2, GP_ds3, Dataset_massflowrate, Dataset_cost, params)
% OBJECTIVEFUN_CEI Computes the Constrained Expected Improvement (CEI).
%
% Maximizes the Expected Improvement (EI) of power consumption, penalized by 
% the Probability of Feasibility (PoF) for surge limits and a Gaussian 
% tracking penalty for the mass flow reference.

    T_row = T(:).';   % Ensure row vector

    % ================= GP Predictions =================
    [mu_J,  std_J]  = gp_predict(GP_Power, T_row);
    [mu_g1, std_g1] = gp_predict(GP_Pdis,  T_row, true); % isP4 = true
    [mu_g2, std_g2] = gp_predict(GP_ds1,   T_row);
    [mu_g3, std_g3] = gp_predict(GP_ds2,   T_row);
    [mu_g4, std_g4] = gp_predict(GP_ds3,   T_row);

    % ================= Expected Improvement (EI) =================
    safe_idxs = Dataset_massflowrate >= massflow_ref;
    if ~any(safe_idxs)
        bestJ = 1e6; % Force exploration if no feasible points exist
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
    % Vectorized for scalability
    mu_dts  = [mu_g2, mu_g3, mu_g4];
    std_dts = max([std_g2, std_g3, std_g4], 1e-12);
    
    % Conservative Lower Confidence Bound (LCB)
    LCB_dts = mu_dts - params.betaSafe .* std_dts;
    
    % Calculate PoF for all constraints simultaneously
    PoF_arr = normcdf_local( (LCB_dts - params.surge_safe_margin) ./ std_dts );
    PoF_arr = max(min(PoF_arr, 1.0), 0.0); % Numerical safety
    
    PoF = prod(PoF_arr); % Combined probability


    % ================= Tracking Penalty =================
        e = (mu_g1 - massflow_ref) / params.cei_track_tol;
        P_track_g1 = 1 / (e.^2 + 1); % The new Scaled Inverse Quadratic
        
    % ================= Final Acquisition =================
    % The GA minimizes the function, so we return the negative of the CEI
    J = - (P_track_g1 .* EI .* PoF);
end

% --- Local Math Helpers ---
function y = normpdf_local(x)
    y = 1./sqrt(2*pi) .* exp(-0.5*x.^2);
end

function y = normcdf_local(x)
    y = 0.5 * (1 + erf(x ./ sqrt(2)));
end
function J = objectiveFun_UCB_series(T, P_dis_ref, GP_Power, GP_Pdis, GP_ds1, GP_ds2, params)

    T_row = T(:).';   % Ensure row vector

    % GP Predictions
    [mu_P,  sigma_P]  = gp_predict(GP_Power, T_row);
    [mu_p,  sigma_g1] = gp_predict(GP_Pdis, T_row, true); % isP4 = true (linear trend)
    [mu_g2, sigma_g2] = gp_predict(GP_ds1, T_row);
    [mu_g3, sigma_g3] = gp_predict(GP_ds2, T_row);

    J_power = (mu_P - params.KappaUCB * sigma_P);

    J_track = (mu_p - P_dis_ref)^2;

    mu_dts    = [mu_g2, mu_g3];
    sigma_dts = [sigma_g2, sigma_g3];
    
    lcb_dts = mu_dts - params.betaSafe * sigma_dts;
    phi = -log(max(lcb_dts - params.surge_safe_margin, 1e-12));
    J_safe = sum(phi);

    J = params.w_power * J_power + params.w_pdis * J_track + params.w_surge * J_safe;
end
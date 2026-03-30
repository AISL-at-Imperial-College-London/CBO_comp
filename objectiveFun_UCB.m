function J = objectiveFun_UCB(T, massflow_ref, GP_Power, GP_Pdis, GP_ds1, GP_ds2, GP_ds3, params)
% OBJECTIVEFUN Computes the composite acquisition function cost.
%
% This function balances power minimization (with uncertainty exploration),
% reference tracking, and log-barrier safety penalties for surge avoidance.

    T_row = T(:).';   % Ensure original torque space is a row vector

    [mu_P,  sigma_P]  = gp_predict(GP_Power, T_row);
    [mu_p,  sigma_g1] = gp_predict(GP_Pdis, T_row, true);
    [mu_g2, sigma_g2] = gp_predict(GP_ds1, T_row);
    [mu_g3, sigma_g3] = gp_predict(GP_ds2, T_row);
    [mu_g4, sigma_g4] = gp_predict(GP_ds3, T_row);

    J_power = (mu_P - params.KappaUCB * sigma_P);

    J_track = (mu_p - massflow_ref)^2;


    mu_dts    = [mu_g2, mu_g3, mu_g4];
    sigma_dts = [sigma_g2, sigma_g3, sigma_g4];
    lcb_dts = mu_dts - params.betaSafe * sigma_dts;
    phi = -log(max(lcb_dts - params.surge_safe_margin, 1e-12));

    J_safe = sum(phi);


    J = params.w_power * J_power + params.w_pdis * J_track + params.w_surge * J_safe;
end
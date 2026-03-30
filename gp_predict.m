function [mu, sigma] = gp_predict(gp, x, isP4)
% GP_PREDICT Predicts mean and std for standard and linear-trend GPs.
    if nargin < 3, isP4 = false; end
    
    x = x(:).';
    x_scaled = (x - gp.lb) ./ gp.span;
    k_star = se_ard_kernel(gp.X, x_scaled, gp.ell, gp.sigma_f);
    
    if isP4
        Phi = [1, x];
        m = Phi * gp.beta;
        mu_scaled = k_star.' * gp.alphaVec;
        v = gp.L \ k_star;
        var_scaled = max(gp.sigma_f^2 - v.'*v, 0);
        
        mu    = m + gp.r_std * mu_scaled;
        sigma = abs(gp.r_std) * sqrt(var_scaled);
    else
        mu_scaled = k_star.' * gp.alphaVec;
        v = gp.L \ k_star;
        var_scaled = max(gp.sigma_f^2 - v.'*v, 0);
        
        mu    = gp.y_mean + gp.y_std * mu_scaled;
        sigma = abs(gp.y_std) * sqrt(var_scaled);
    end
end
function gp = make_gp_p4(X_scaled, X_orig, y, lb, span)
% MAKE_GP_P4 GP with linear trend mean and SE-ARD kernel on residuals.
    [N, d] = size(X_scaled);
    
% Fit linear trend dynamically for any number of compressors
    Phi = [ones(N,1), X_orig];
    beta = Phi \ y;
    m_train = Phi * beta;      
    r       = y - m_train;     % residuals
    
    r_std = std(r);
    if r_std == 0, r_std = 1; end
    r_scaled = r / r_std;
    
    ell     = 0.25 * ones(d,1); 
    sigma_f = 1.0;
    alpha   = 1e-4;
    
    K = se_ard_kernel(X_scaled, X_scaled, ell, sigma_f);
    K = K + alpha * eye(N);
    L = chol(K, 'lower');
    alpha_vec = L'\(L\r_scaled);
    
    gp.X        = X_scaled;
    gp.beta     = beta;
    gp.r_std    = r_std;
    gp.ell      = ell;
    gp.sigma_f  = sigma_f;
    gp.alpha    = alpha;
    gp.L        = L;
    gp.alphaVec = alpha_vec;
    gp.lb       = lb(:).';
    gp.span     = span(:).';
end
function gp = make_gp(X_scaled, y, lb, span, kernelType)
% MAKE_GP Constructs a standard GP surrogate model.
    [~, d] = size(X_scaled);
    y_mean = mean(y);
    y_std  = std(y);
    if y_std == 0, y_std = 1; end
    
    y_scaled = (y - y_mean)/y_std;
    ell     = 0.5 * ones(d,1);
    sigma_f = 1.0;
    alpha   = 1e-4;
    
    K = se_ard_kernel(X_scaled, X_scaled, ell, sigma_f);
    K = K + alpha*eye(size(X_scaled,1));
    L = chol(K, 'lower');
    alpha_vec = L'\(L\y_scaled);
    
    gp.kernelType = kernelType;
    gp.X        = X_scaled;
    gp.y_mean   = y_mean;
    gp.y_std    = y_std;
    gp.ell      = ell;
    gp.sigma_f  = sigma_f;
    gp.alpha    = alpha;
    gp.L        = L;
    gp.alphaVec = alpha_vec;
    gp.lb   = lb(:).';
    gp.span = span(:).';
end
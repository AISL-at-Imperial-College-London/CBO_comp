function K = se_ard_kernel(X1, X2, ell, sigma_f)
% SE_ARD_KERNEL Squared Exponential Kernel with Automatic Relevance Determination.
    X1n = X1 ./ ell.';
    X2n = X2 ./ ell.';
    X1sq = sum(X1n.^2, 2);
    X2sq = sum(X2n.^2, 2).';
    D2   = X1sq + X2sq - 2*(X1n*X2n');
    K = (sigma_f^2)*exp(-0.5*D2);
end
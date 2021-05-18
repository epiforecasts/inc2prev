  covariates = update_gp(PHI, M, L, alpha, rho, eta, 0);
  covariates = cumulative_sum(covariates);


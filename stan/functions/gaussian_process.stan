// all from here: 
// https://github.com/epiforecasts/EpiNow2/blob/master/inst/stan/functions/gaussian_process.stan
// eigenvalues for approximate hilbert space gp
// see here for details: https://arxiv.org/pdf/2004.11408.pdf
real lambda(real L, int m) {
  real lam;
  lam = ((m*pi())/(2*L))^2;
  return lam;
}
// eigenfunction for approximate hilbert space gp
// see here for details: https://arxiv.org/pdf/2004.11408.pdf
vector phi(real L, int m, vector x) {
  vector[rows(x)] fi;
  fi = 1/sqrt(L) * sin(m*pi()/(2*L) * (x+L));
  return fi;
}
// spectral density of the exponential quadratic kernal
real spd_se(real alpha, real rho, real w) {
  real S;
  S = (alpha^2) * sqrt(2*pi()) * rho * exp(-0.5*(rho^2)*(w^2));
  return S;
}
// spectral density of the Matern 3/2 kernel
real spd_matern(real alpha, real rho, real w) {
  real S;
  S = 4*alpha^2 * (sqrt(3)/rho)^3 * 1/((sqrt(3)/rho)^2 + w^2)^2;
  return S;
}
// setup approximate gaussian process
matrix setup_gp(int M, real L, int dimension) {
  vector[dimension] time;
  matrix[dimension, M] PHI;
  real half_dim = dimension / 2.0;
  for (s in 1:dimension) {
    time[s] = (s - half_dim) / half_dim;
  }
  for (m in 1:M){
    PHI[,m] = phi(L, m, time);
  }
  return(PHI);
}
// update gaussian process using spectral densities
vector update_gp(matrix PHI, int M, real L, real alpha,
                 real rho, vector eta, int type) {
  vector[M] diagSPD;    // spectral density
  vector[M] SPD_eta;    // spectral density * noise
  int noise_terms = rows(PHI);
  vector[noise_terms] noise = rep_vector(1e-8, noise_terms);
  real unit_rho = rho / noise_terms;
  // GP in noise - spectral densities
  if (type == 0) {
    for(m in 1:M){
      diagSPD[m] = sqrt(spd_se(alpha, unit_rho, sqrt(lambda(L, m))));
    }
  }else if (type == 1) {
    for(m in 1:M){
      diagSPD[m] = sqrt(spd_matern(alpha, unit_rho, sqrt(lambda(L, m))));
    }
  }
  SPD_eta = diagSPD .* eta;
  noise = noise + PHI[,] * SPD_eta;
  return(noise);
}

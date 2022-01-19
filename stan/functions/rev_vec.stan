vector rev_vec(vector vec) {
  int l = num_elements(vec);
  vector[l] rvec;
  for (i in 1:l) {
    rvec[i] = vec[l - i + 1];
  }
  return(rvec);
}

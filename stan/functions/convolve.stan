vector convolve(vector vec1, vector vec2) {
  int lvec1 = num_elements(vec1);
  int lvec2 = num_elements(vec2);
  vector[lvec1] conv = rep_vector(1e-8, lvec1);
  for (s in 1:lvec1) {
    conv[s] += dot_product(vec1[max(1, (s - lvec2 + 1)):s], 
                           tail(vec2, min(lvec2, s)));
  }
  return(conv);
}


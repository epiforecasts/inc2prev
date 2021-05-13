#' Posterior samples for the probability of being detected for COVID-19
#'
#' @description An example data set indicating the format required for
#' samples of the probability of being detected. Each row should represent a
#' sample with the first variable representing the sample id and the second
#' variable should indicate the time at which the sample applies. The following
#' variables should contain the probability of detection on each day starting
#' at day 0 and  until the limit of detection is reached with no missing days.
#'
#' @details These estimates were extracted from: "Estimating the effectiveness
#'  of routine asymptomatic PCR testing at different frequencies for the
#'  detection of SARS-CoV-2 infections" by Hellewell et al.
#'
#' If using these estimates in a publication please cite:
#
#' Hellewell, J., Russell, T.W., The SAFER Investigators and Field Study Team
#' et al. Estimating the effectiveness of routine asymptomatic PCR testing at
#' different frequencies for the detection of SARS-CoV-2 infections. BMC Med
#' 19, 106 (2021). https://doi.org/10.1186/s12916-021-01982-x
#
#' @return A data frame with samples of the probability of detection in the
#' appropriate format.
"covid19_prob_detectable"

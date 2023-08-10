my_cols <- c("GLM-EIV (accelerated)" = "firebrick3",
             "Thresholding (w/ oracle threshold)" = "dodgerblue3",
             "Unimodal mixture" = "orchid4", "GLM-EIV" = "dimgrey")

transform_metric_df <- function(metric_df) {
  out <- metric_df |>
    dplyr::filter(metric %in% c("bias", "ci_width", "coverage", "mse", "time")) |>
    mutate(metric_fct = factor(metric, levels = c("bias", "mse", "coverage", "ci_width", "time"),
                               labels = c("Bias", "MSE", "Coverage", "CI width", "Time (s)")),
           Method = factor(method, levels = c("glmeiv_slow", "glmeiv_fast", "thresholding", "unimodal_mixture"),
                           labels = c("GLM-EIV", "GLM-EIV (accelerated)", "Thresholding (w/ oracle threshold)", "Unimodal mixture"))) |>
    arrange(Method)
  if ("fam_str" %in% colnames(out)) {
    out <- out |> mutate(fam_str = factor(fam_str, levels = c("nb_theta_unknown", "nb_theta_known", "poisson"),
                         labels = c("NB (theta est.)", "NB (theta known)", "Poisson")))
    
  }
  return(out)
}
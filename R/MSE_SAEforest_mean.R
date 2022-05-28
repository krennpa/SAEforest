# MSE bootstrap for means under unit-level covariate data ---------------------------------

MSE_SAEforest_mean <- function(Y,
                               X,
                               dName,
                               smp_data,
                               mod,
                               ADJsd,
                               pop_data,
                               B = 100,
                               initialRandomEffects = 0,
                               ErrorTolerance = 0.0001,
                               MaxIterations = 25,
                               ...) {

  rand_struc <- paste0(paste0("(1|", dName), ")")
  domains <- t(unique(pop_data[dName]))
  in_samp <- domains %in% t(unique(smp_data[dName]))
  N_i <- as.numeric(table(pop_data[[dName]]))

  # Prepare data for sampling
  ran_obj <- ran_comp(Y = Y, smp_data = smp_data, mod = mod, ADJsd = ADJsd, dName = dName)
  ran_effs <- ran_obj$ran_effs
  forest_res <- ran_obj$forest_res
  smp_data <- ran_obj$smp_data

  pred_vals <- predict(mod$Forest, pop_data)$predictions
  pred_mat <- matrix(pred_vals, nrow = length(pred_vals), ncol = B)

  block_sample_e <- function(x) {
    block_sample(
      domains = domains, in_samp = in_samp,
      smp_data = smp_data, dName = dName,
      pop_data = pop_data, forest_res = forest_res
    )
  }

  sample_ui <- function(x) {
    rep(sample(ran_effs,
      size = length(N_i),
      replace = TRUE
    ), N_i)
  }

  e_ij <- matrix(NA, nrow = length(pred_vals), ncol = B)
  e_ij <- apply(e_ij, 2, block_sample_e)

  u_i <- apply(pred_mat, 2, sample_ui)

  smp_data$forest_res <- NULL

  # combine to y_star
  y_star <- pred_mat + u_i + e_ij

  indi_agg <- rep(1:length(N_i), N_i)
  my_agg <- function(x) {
    tapply(x, indi_agg, mean)
  }
  tau_star <- apply(y_star, my_agg, MARGIN = 2)

  # get bootstrap samples
  boots_sample <- vector(mode = "list", length = B)

  for (i in 1:B) {
    pop_data$y_star <- y_star[, i]
    boots_sample[[i]] <- sample_select(pop_data, smp = smp_data, dName = dName)
  }

  # use bootstrap samples to estimate
  my_estim_f <- function(x) {
    point_mean(
      Y = x$y_star, X = x[, colnames(X)], dName = dName, smp_data = x,
      pop_data = pop_data, initialRandomEffects = initialRandomEffects,
      ErrorTolerance = ErrorTolerance, MaxIterations = MaxIterations, ...
    )[[1]][, 2]
  }

  tau_b <- pbapply::pbsapply(boots_sample, my_estim_f)

  MSE_estimates <- rowMeans((tau_star - tau_b)^2)

  MSE_estimates <- data.frame(unique(pop_data[dName]), Mean = MSE_estimates)
  rownames(MSE_estimates) <- NULL

  return(MSE_estimates)
}

# MSE bootstrap for means under aggregated covariate data ---------------------------------
MSE_SAEforest_aggOOB <- function(Y,
                                 X,
                                 dName,
                                 smp_data,
                                 mod,
                                 ADJsd,
                                 Xpop_agg,
                                 B = 100,
                                 popnsize,
                                 initialRandomEffects,
                                 ErrorTolerance,
                                 MaxIterations,
                                 ...) {

  dom <- smp_data[[dName]]
  in_dom <- as.character(unique(smp_data[[dName]]))
  total_dom <- as.character(mod$Indicators[[dName]])
  p <- dim(X)[2]

  Freq_n <- data.frame(table(dom))
  colnames(Freq_n)[1] <- dName
  n_i <- Freq_n$Freq
  freq <- dplyr::left_join(popnsize, Freq_n, by = dName)

  freq[, 3][is.na(freq[, 3])] <- 0
  freq <- freq[match(total_dom, popnsize[[dName]]), ]

  rd <- freq[, 2] - freq[, 3]
  N_i <- freq[, 2]

  # Prepare data for sampling
  ran_obj <- ran_comp(
    Y = Y, smp_data = smp_data, mod = mod$MERFmodel, ADJsd = ADJsd,
    dName = dName
  )
  ran_effs <- ran_obj$ran_effs
  forest_res <- ran_obj$forest_res
  smp_data <- ran_obj$smp_data

  pred_mat <- matrix(mod$MERFmodel$Forest$predictions, nrow = length(dom), ncol = B)

  random_eff <- tapply(mod$ModifiedSet$u_ij, mod$ModifiedSet[[dName]], mean)
  mu_t <- mod$Indicators$Mean - random_eff

  mu_pred <- matrix(mu_t, nrow = length(total_dom), ncol = B)

  e_ij <- matrix(sample(forest_res, size = length(pred_mat), replace = TRUE),
    nrow = length(dom), ncol = B, byrow = FALSE
  )

  e_i_mean <- as.matrix(aggregate(. ~ dom, data.frame(dom, e_ij), FUN = mean)[, -1])

  u_i <- apply(mu_pred, 2, function(x) {
    sample(ran_effs, size = length(total_dom), replace = TRUE)
  })

  u_ij <- apply(u_i[total_dom %in% in_dom, ], 2, function(x) {
    rep(x, n_i)
  })

  # combine to y_star
  y_star <- pred_mat + u_ij + e_ij

  samp_erd <- function(x) {
    sample((forest_res / ADJsd) * sqrt(ADJsd^2 / x), size = B, replace = TRUE)
  }
  erd_mean <- t(sapply(rd, samp_erd))

  insamp_ei <- e_i_mean * n_i / N_i[total_dom %in% in_dom] +
    erd_mean[total_dom %in% in_dom, ] * (rd / N_i)[total_dom %in% in_dom]

  tau_star <- mu_pred + u_i + erd_mean
  tau_star[total_dom %in% in_dom, ] <- mu_pred[total_dom %in% in_dom, ] +
    u_i[total_dom %in% in_dom, ] + insamp_ei

  # use bootstrap samples to estimate
  my_estim_f <- function(x) {
    point_meanAGG(
      Y = x, X = X, dName = dName, smp_data = smp_data,
      Xpop_agg = Xpop_agg, initialRandomEffects = initialRandomEffects,
      ErrorTolerance = ErrorTolerance, MaxIterations = MaxIterations,
      OOsample_obs = mod$OOsample_obs, ADDsamp_obs = mod$ADDsamp_obs,
      w_min = mod$w_min, wSet = mod$wSet, verbose = FALSE, ...
    )[[1]]$Mean
  }

  tau_b <- pbapply::pbapply(y_star, 2, my_estim_f)

  MSE_estimates <- rowMeans((tau_star - tau_b)^2)

  MSE_estimates <- data.frame(mod$Indicators[dName], Mean = MSE_estimates)
  rownames(MSE_estimates) <- NULL

  return(MSE_estimates)
}

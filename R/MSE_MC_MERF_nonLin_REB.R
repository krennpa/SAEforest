MSE_MC_MERF_nonLin_REB <- function(Y, X, dName, threshold, smp_data, mod, ADJsd, pop_data, B=100, B_point,
                                     initialRandomEffects = 0, ErrorTolerance = 0.0001,
                                     MaxIterations = 25, custom_indicator,...){

  forest_m1 <- mod
  rand_struc = paste0(paste0("(1|",dName),")")
  domains <- t(unique(pop_data[dName]))

  n_i <- as.numeric(table(pop_data[[dName]]))

  boots_pop <- vector(mode="list",length = B)
  boots_pop <- sapply(boots_pop,function(x){pop_data},simplify =FALSE)

  # DATA PREP ________________________________________________
  ran_obj <- ran_comp(Y=Y, smp_data = smp_data, mod=forest_m1, ADJsd = ADJsd, dName = dName)
  ran_effs <- ran_obj$ran_effs
  forest_res <- ran_obj$forest_res
  smp_data <- ran_obj$smp_data

  pred_vals<- predict(forest_m1$Forest, pop_data)$predictions
  my_pred_f <- function(x){pred_vals}
  pred_t <- vector(mode="list", length = B)
  pred_t <- sapply(pred_t, my_pred_f,simplify = FALSE)


  # Errors
  samp_e <- function(x){sample(forest_res, size=sum(n_i), replace = TRUE)}
  e_ij <- sapply(boots_pop, samp_e, simplify = FALSE)

  smp_data$forest_res <- NULL

  u_i <- replicate(length(boots_pop),rep(sample(ran_effs, size = length(t(unique(pop_data[dName]))),
                                                replace=TRUE), n_i), simplify = FALSE)

  # combine
  y_star_u_star <-  Map(function(x,y,z){x+y+z}, pred_t, e_ij, u_i)

  if(is.numeric(threshold)){
    thresh <- sapply(y_star_u_star, function(x){threshold}, simplify = FALSE)
  }
  if(is.null(threshold)){
    thresh <- sapply(y_star_u_star, function(x){0.6*median(x, na.rm=TRUE)}, simplify = FALSE)
  }
  if(is.function(threshold)){
    thresh <- sapply(y_star_u_star, threshold, simplify = FALSE)
  }

  boots_pop <- Map(cbind, boots_pop, "y_star_u_star" = y_star_u_star, "thresh"=thresh)


  my_agg <-  function(x){tapply(x[["y_star_u_star"]], x[[dName]], calc_indicat, threshold = unique(x[["thresh"]]), custom = custom_indicator)}
  tau_star <- sapply(boots_pop, my_agg, simplify = FALSE)
  comb <- function(x){matrix(unlist(x), nrow = length(n_i), byrow = TRUE)}
  tau_star <- sapply(tau_star, comb, simplify = FALSE)


  # THINK ABOUT SEED
  sample_b <- function(x){sample_select(x,smp=smp_data,times=1, dName = dName)}
  boots_sample <- sapply(boots_pop, sample_b)


  # USE BOOTSTRAP SAMPLE TO ESITMATE

  my_estim_f <- function(x){point_MC_nonLin(Y = x$y_star_u_star, X = x[,colnames(X)], dName = dName, threshold = threshold, smp_data = x, pop_data = pop_data,
                                         initialRandomEffects = initialRandomEffects, ErrorTolerance = ErrorTolerance, B_point = B_point,
                                         MaxIterations = MaxIterations, custom_indicator=custom_indicator,...)[[1]][,-1]}

  tau_b <- pbapply::pbsapply(boots_sample, my_estim_f,simplify = FALSE)

  mean_square <- function(x,y){(x-y)^2}

  Mean_square_B <- mapply(mean_square, tau_b,tau_star, SIMPLIFY = FALSE)

  MSE_estimates <- Reduce('+',Mean_square_B)/length(Mean_square_B)
  MSE_estimates_out <- data.frame(unique(pop_data[dName]), MSE_estimates)
  colnames(MSE_estimates_out) <- c(dName, colnames(MSE_estimates))
  rownames(MSE_estimates_out) <- NULL

  #__________________________
  return(MSE_estimates_out)

}

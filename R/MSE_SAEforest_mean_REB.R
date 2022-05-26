MSE_SAEforest_mean_REB <- function(Y, X, dName, smp_data, mod, ADJsd, pop_data, B=100,
                                   initialRandomEffects = 0, ErrorTolerance = 0.0001,
                                   MaxIterations = 25, ...){

  forest_m1 <- mod
  rand_struc = paste0(paste0("(1|",dName),")")
  domains <- t(unique(pop_data[dName]))
  boots_pop <- vector(mode="list",length = B)
  boots_pop <- sapply(boots_pop,function(x){pop_data},simplify =FALSE)
  n_i <- as.numeric(table(pop_data[[dName]]))

 # DATA PREP  ________________________________________________
  ran_obj <- ran_comp(Y=Y, smp_data = smp_data, mod=forest_m1, ADJsd = ADJsd, dName = dName)
  ran_effs <- ran_obj$ran_effs
  forest_res <- ran_obj$forest_res
  smp_data <- ran_obj$smp_data

  pred_vals<- predict(forest_m1$Forest, pop_data)$predictions
  my_pred_f <- function(x){pred_vals}
  pred_t <- vector(mode="list", length = B)
  pred_t <- sapply(pred_t, my_pred_f,simplify = FALSE)

  # Errors
  block_sample <- function(x){
    in_samp <- domains %in% t(unique(smp_data[dName]))
    block_err <- vector(mode="list",length = length(domains))

    for (idd in which(in_samp)){
      block_err[[idd]] <- sample(forest_res[smp_data[dName]==domains[idd]],size = sum(pop_data[dName] == domains[idd]),replace=TRUE)
    }
    if (sum(in_samp)!= length(domains)){
      for (idd in which(!in_samp)){
        block_err[[idd]] <- sample(forest_res, size = sum(pop_data[dName] == domains[idd]), replace=TRUE)
      }
    }
    return(unlist(block_err))
  }

  e_ij <- sapply(boots_pop, block_sample, simplify = FALSE)

  smp_data$forest_res <- NULL

  u_i <- replicate(length(boots_pop),rep(sample(ran_effs, size = length(n_i),
                                                replace=TRUE), n_i), simplify = FALSE)

  # combine
  y_star_u_star <-  Map(function(x,y,z){x+y+z}, pred_t, e_ij, u_i)

  boots_pop <- Map(cbind, boots_pop, "y_star_u_star" = y_star_u_star)

  agg_form <- formula(paste("y_star_u_star ~ ", paste0(dName)))
  my_agg <- function(x){aggregate(agg_form, x,  mean)[,2]}
  tau_star <- sapply(boots_pop, my_agg)

  # THINK ABOUT SEED
  sample_b <- function(x){sample_select(x, smp = smp_data, times = 1, dName = dName)}
  boots_sample <- sapply(boots_pop, sample_b)


  # USE BOOTSTRAP SAMPLE TO ESITMATES
  my_estim_f <- function(x){point_mean(Y = x$y_star_u_star, X = x[,colnames(X)], dName = dName, smp_data = x,
                                       pop_data = pop_data, initialRandomEffects = initialRandomEffects,
                                       ErrorTolerance = ErrorTolerance, MaxIterations = MaxIterations)[[1]][,2]}

  tau_b <- pbapply::pbsapply(boots_sample, my_estim_f)

  MSE_estimates <- rowMeans((tau_star - tau_b)^2)

  MSE_estimates <- data.frame(unique(pop_data[dName]), Mean=MSE_estimates)
  rownames(MSE_estimates) <- NULL

  #___________________________

  return(MSE_estimates)

}



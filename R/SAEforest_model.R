

SAEforest_model <- function(Y, X, dName, smp_data, pop_data, MSE = "none", meanOnly = TRUE,
                            aggData =FALSE, smearing = TRUE, popnsize = NULL,
                            importance = "impurity", OOsample_obs = 25, ADDsamp_obs=0,
                            w_min=3, B=100, B_adj = 100, B_MC=100, threshold = NULL,
                            custom_indicator =NULL, initialRandomEffects = 0,
                            ErrorTolerance = 0.0001, MaxIterations = 25, na.rm =TRUE,...){

  input_checks_model(Y =  Y, X =X, dName = dName, smp_data = smp_data, pop_data = pop_data,
                     MSE = MSE, meanOnly = meanOnly, aggData =aggData, smearing = smearing,
                     popnsize = popnsize, importance = importance, OOsample_obs = OOsample_obs,
                     ADDsamp_obs=ADDsamp_obs, w_min=w_min, B=B, B_adj = B_adj, B_MC=B_MC,
                     threshold = threshold, custom_indicator =custom_indicator,
                     initialRandomEffects = initialRandomEffects, ErrorTolerance = ErrorTolerance,
                     MaxIterations = MaxIterations, na.rm =na.rm)

  if(meanOnly == TRUE || aggData == TRUE){
      return_obj <- SAEforest_mean(Y = Y, X = X, dName = dName, smp_data = smp_data, pop_data =pop_data,
                       MSE = MSE, aggData =aggData, popnsize = popnsize, OOsample_obs = OOsample_obs,
                       ADDsamp_obs=ADDsamp_obs, w_min=w_min, importance = importance,
                       initialRandomEffects = initialRandomEffects, ErrorTolerance = ErrorTolerance,
                       MaxIterations = MaxIterations,  B=B, B_adj = B_adj, na.rm =na.rm,...)

      return(return_obj)
  }

  if(meanOnly == FALSE){
      return_obj <- SAEforest_nonLin(Y = Y, X = X, dName = dName, smp_data = smp_data, pop_data = pop_data,
                       smearing = smearing, MSE = MSE, importance =importance,
                       initialRandomEffects = initialRandomEffects, ErrorTolerance = ErrorTolerance,
                       MaxIterations = MaxIterations, B=B, B_adj=B_adj, B_MC=B_MC, threshold = threshold,
                       custom_indicator =custom_indicator, na.rm = na.rm, ...)

  return(return_obj)
  }
}


# General Support Functions ----

# Recommeded to add by the check function
#' @importFrom stats aggregate as.formula complete.cases dist formula median na.omit optimize predict quantile rnorm sd


calc_indicat <- function(Y, threshold, custom){

  hcr_function <- function(y,threshold){
    mean(y < threshold, na.rm = TRUE)
  }
  qsr_function <- function(y){
    sum(y[(y > quantile(y,0.8, na.rm = TRUE))]) / sum(y[(y < quantile(y,0.2, na.rm = TRUE))])
  }
  pgap_function <- function(y,threshold) {
    y[y<0] <- NA
    mean((y < threshold)*(threshold - y) / threshold, na.rm = TRUE)
  }
  gini_function <- function(y){
    y[y<0] <- NA
    ineq::Gini(y, na.rm = TRUE)
  }
  quant_preds <- quantile(Y, prob=c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)
  mean_est <- mean(Y, na.rm = TRUE)
  Gini_est <- gini_function(y = Y)
  Hcr_est <- hcr_function(y = Y, threshold = threshold)
  Qsr_est <- qsr_function(y = Y)
  Pgap_est <- pgap_function(y = Y, threshold = threshold)

  indicators <- cbind(mean_est, t(quant_preds), Gini_est, Hcr_est, Pgap_est, Qsr_est)

  colnames(indicators) <- c("Mean","Quant10","Quant25","Median","Quant75",
                          "Quant90","Gini","Hcr","Pgap","Qsr")

  if(!is.null(custom)){
    custom_ind <- unlist(lapply(custom, function(f) f(Y, threshold)))
    indicators <- cbind(indicators,t(custom_ind))
  }

  return(indicators)

}

sae_specs <- function(dName,cns,smp){
  in_dom<- unique(smp[[dName]])
  total_dom <- unique(cns[[dName]])
  OOsamp <- !total_dom %in% in_dom

  return(list(
    n_surv = length(smp[[dName]]),
    n_pop = length(cns[[dName]]),
    n_out = sum(OOsamp),
    n_in = length(in_dom),
    n_total = length(total_dom),
    n_smp = table(smp[[dName]]),
    n_cns = table(cns[[dName]])))
}


# DEFINE AND USE SAMPLE SELECT WRAPPER FUNCTION in MSEs ----
sample_select <- function(pop, smp, dName){
  # pop.............. the population or census data
  # smp.............. the sample data
  # set_seed......... set seed for reproduceability

  smpSizes <- table(smp[dName])
  smpSizes <- data.frame(smpidD = as.character(names(smpSizes)), n_smp = as.numeric(smpSizes),
                         stringsAsFactors = FALSE)

  smpSizes <- dplyr::left_join(data.frame(idD = as.character(unique(pop[[dName]]))),
                               smpSizes, by = c("idD" = "smpidD"))

  smpSizes$n_smp[is.na(smpSizes$n_smp)] <- 0

  splitPop <- split(pop, pop[[dName]])

  stratSamp <- function(dfList, ns) {
    do.call(rbind, mapply(dfList, ns, FUN = function(df, n) {
      popInd <- seq_len(nrow(df))
      sel <- base::sample(popInd, n, replace = FALSE)
      df[sel, ]
    }, SIMPLIFY = F))
  }

  samples <- stratSamp(splitPop, smpSizes$n_smp)

  rm(splitPop)
  return(samples)
}

ran_comp <- function(mod, smp_data, Y, dName, ADJsd){

  forest_res1 <- Y - predict(mod$Forest, smp_data)$predictions
  smp_data$forest_res <- forest_res1

  # Random Effects
  formRF <- formula(paste("forest_res ~", paste0(dName)))
  ran_effs1 <- aggregate(data=smp_data, formRF, FUN=mean)
  colnames(ran_effs1) <- c(dName,"r_bar")

  smp_data <- dplyr::left_join(smp_data,ran_effs1,by = dName)
  smp_data$forest_eij <- smp_data$forest_res-smp_data$r_bar

  # prepare for sampling
  forest_res <- smp_data$forest_eij
  forest_res<-(forest_res/sd(forest_res))*ADJsd

  # CENTER
  forest_res <- forest_res-mean(forest_res)

  # prepare for sampling
  ran_effs <- ran_effs1$r_bar
  ran_effs <- (ran_effs/sd(ran_effs))*mod$RanEffSD

  # CENTER
  ran_effs <- ran_effs-mean(ran_effs)

  return(list(forest_res = forest_res,
              ran_effs = ran_effs,
              smp_data = smp_data))
}

ran_comp_wild <- function(mod, smp_data, Y, dName, ADJsd){

  forest_res <- Y - mod$Forest$predictions-predict(mod$EffectModel, smp_data)
  forest_res <- forest_res-mean(forest_res)

  # Random Effects
  ran_effs <- unique(predict(mod$EffectModel, smp_data))
  ran_effs <- (ran_effs/sd(ran_effs))* mod$RanEffSD
  ran_effs <- ran_effs-mean(ran_effs)

  return(list(forest_res = forest_res,
              ran_effs = ran_effs))
}

# Errors
block_sample <- function(domains, in_samp, smp_data, dName, pop_data, forest_res){

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


wild_errors <- function(x, mod, smp_data, forest_res){
  fitted_s <- predict(mod$Forest, smp_data)$predictions + predict(mod$EffectModel, smp_data)
  indexer <- vapply(x, function(x) {which.min(abs(x - fitted_s))},
                    FUN.VALUE = integer(1))
  # superpopulation individual errors
  eps <- forest_res[indexer]
  wu <- sample(c(-1,1),size = length(eps), replace = TRUE)
  eps <- abs(eps) * wu

  return(eps)}


get_thresh <- function(x,threshold){
  if(is.numeric(threshold)){
    thresh <- threshold
  }
  if(is.null(threshold)){
    thresh <- 0.6*median(x, na.rm=TRUE)
  }
  if(is.function(threshold)){
    thresh <- threshold(x)
  }
  return(thresh)
}


# EMPIRICAL LIKELIHOOD SUPPORT FUNS ----
#
elm <- function(x, mu = rep(0,ncol(x)), lam = t(rep(0,ncol(x))), maxit = 25,
                gradtol = 1e-7, svdtol = 1e-9, itertrace = 0){

  n = nrow(x)
  p = ncol(x)

  if(ncol(x) != length(mu)){stop('Need size(mu) = size(x(i,:)')}

  if(length(lam) != length(mu)){stop('Lam must be sized like mu')}

  if(gradtol <  1e-16){gradtol = 1e-16}
  if(svdtol  <  1e-8){gradtol   = 1e-8}

  z = x - t(matrix(mu, nrow = p, ncol= n))

  # Step weights: inner search with Newton step, then if necessary gets smaller
  # and more parallel to the gradient

  newton_wts     =   c(1/3 ^c(0:3), rep(0,12))
  gradient_wts   =   .5 ^c(0:15)
  gradient_wts   =   (gradient_wts^2 - newton_wts^2)^0.5

  gradient_wts[12:16] = gradient_wts[12:16] / (.1^-(1:5))

  #  Outer loop of iteration. Newton iteration

  nits  = 0
  gsize = gradtol + 1.0

  while(nits < maxit & gsize > gradtol){

    arg  = 1 + z %*% lam
    wts1 = plog(arg ,1/n,1 )
    wts2 = (- plog( arg,1/n,2 ))^0.5

    grad = -z * matrix(wts1, nrow = nrow(wts1), ncol = p)
    grad = colSums(grad)
    gsize= mean(abs(grad))

    hess = z * matrix(wts2, nrow = nrow(wts2), ncol = p) # matrix sqrt of hessian

    hu  = svd(hess)$u
    hs  = diag(length(svd(hess)$d)) * svd(hess)$d
    hv  = svd(hess)$v
    dhs = svd(hess)$d

    if(min(dhs) < max(dhs)*svdtol + 1e-128){
      dhs = dhs + svdtol*max(dhs) + 1e-128
    }

    dhs= 1/dhs
    hs = diag(dhs)

    nstep = hv %*% hs %*% t(hu) %*% (wts1/wts2)
    gstep = -grad
    if(sum(nstep^2) < sum(gstep^2)){
      gstep = gstep*sum(nstep^2)^0.5/sum(gstep^2)^0.5
    }

    ologelr = -sum(plog(arg, 1/n))
    ninner  = 0

    for(i in 1:length(newton_wts)){
      nlam    = lam + newton_wts[i] * nstep + gradient_wts[i] * gstep
      nlogelr = -sum(plog(1 + z %*% nlam, 1/n))

      if(nlogelr < ologelr){
        lam    = nlam
        ninner = i
        break
      }
    }

    nits = nits+1
    if(ninner == 0){
      nits = maxit
    }

    if(itertrace == 1){
      print(list(lam,nlogelr,gsize,ninner,nits))
    }
  }

  logelr = nlogelr
  lambda = lam
  hess   = t(hess) %*% hess
  wts    = wts1

  Ergebnisse <- list(logelr, lambda, grad, hess, wts, nits)
  names(Ergebnisse) <- c("logelr", "lambda", "grad", "hess", "wts", "nits")

  return(Ergebnisse)
}


plog <- function(z,eps,d = 0){

  # Pseudo logarithm
  # The pseudo logarithm agrees with log for arguments larger than eps.

  zsize = dim(z)
  out   = c(z)
  low   = out < eps

  if(d == 0){
    out[!low] = log(out[!low])
    out[ low] = log(eps) - 1.5 + 2*out[low]/eps - 0.5*(out[low]/eps)^2
  }else if(d == 1){
    out[!low] = 1/out[!low]
    out[ low] = 2/eps - out[low]/eps^2
  }else if(d == 2){
    out[!low] = -1/out[!low]^2
    out[ low] = -1/eps^2
  }else{
    stop('Unknown option d for pseudologarithm')
  }
  llogz = matrix(out, zsize[1], zsize[2])
  return(llogz)
}

# Wrapper Function

elm_wrapper <- function(X_input_elm,mu_input_elm){
  n = nrow(X_input_elm)
  p = ncol(X_input_elm)

  elm_results <- elm(x = X_input_elm , mu =  mu_input_elm, lam = rep(0, p), maxit = 25, gradtol = 1e-7,
                     svdtol = 1e-9, itertrace = 0)

  logelr = elm_results$logelr
  lambda = elm_results$lambda

  zz    = X_input_elm - t(matrix(mu_input_elm, nrow = p, ncol = n))
  tmp3   = n * (zz %*% lambda + 1)
  # Gewicht mit Nebenbedingungen bestimmen
  prob   = 1/tmp3

  return(list(prob = prob, experWeight = elm_results$wts))
}


sortAlpha <- function(x, dName){
  orderIndicat <- order(as.character(x[[dName]]))
  return_val <- x[orderIndicat,]
  rownames(return_val) <- NULL
  return(return_val)
}

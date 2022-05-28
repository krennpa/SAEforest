# Auxiliary function for MSEs to perform bias correction ---------------------------------

adjust_ErrorSD <- function(Y, X, smp_data, mod, B = 100, ...) {
  pred_OOB <- matrix(mod$Forest$predictions, ncol = B,
                     nrow = length(mod$Forest$predictions), byrow = FALSE)

  e_ij <- Y - predict(mod$Forest, smp_data)$predictions
  e_ij <- e_ij - mean(e_ij)

  y_star_OOB <- pred_OOB + sample(e_ij, size = length(pred_OOB), replace = TRUE)

  my_estim_f2 <- function(x) {
    ranger::ranger(y = x, x = X, data = smp_data, ...)
  }
  my_f_n2 <- pbapply::pbapply(y_star_OOB, 2, my_estim_f2)

  my_pred_f <- function(x) {
    x$predictions
  }
  pred_OOB_star <- sapply(my_f_n2, my_pred_f)

  Adjustment <- (pred_OOB_star - pred_OOB)^2

  outvar <- sqrt(mod$ErrorSD^2 - mean(rowMeans(Adjustment)))

  return(outvar)
}

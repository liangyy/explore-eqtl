# x: n-by-p matrix
# y: n-by-1 vector
# z: n-by-k matrix
fineMapping <- function(x, y, z) {
  y <- .lm.fit(z, y)$residuals
  obj <- susieR::susie(x, y, L = 5, estimate_residual_variance = TRUE, prior_variance = 1, intercept = TRUE, tol = 1e-3)
  return(obj)
}

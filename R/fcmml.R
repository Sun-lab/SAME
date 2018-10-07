#' Fisher score method for estimating the regression coeffcients in a mixture of logisitc regression.
#'
#' @param Y The response variable. Must be binary.
#' @param X The design matrix. Intercept included.
#' @param alpha The initial values of the coefficients for the confounders.
#' @param beta The initial value of the coefficient for the true somatic mutation value.
#' @param T0S The weights vector for the mixture when the true somatic muation equals 0.
#' @param T1S The weights vector for the mixture when the true somatic muation equals 1.
#' @param converged The tolerance for the convergence. Default is 1e-6.
#' @param maxIT The maximal number of the EM iteration times. Default is 200.
#'
#' @return The estimated regression coefficients.









fcmml <- function(Y, X, alpha, beta, T0S, T1S, converged = 1e-6,
                  maxIT = 200){

  n <- nrow(X)
  d <- ncol(X)

  theta_old <- theta <- c(alpha, beta)

  for(i in 1:maxIT){
    alpha <- theta_old[1:d]
    beta <- theta_old[d + 1]

    # weights
    w1_alpha <- exp(X %*% alpha ) / (1 + exp(X %*% alpha ))
    w1_beta <- exp(X %*% alpha + beta) / (1 + exp(X %*% alpha + beta))

    w2_alpha <- exp(X %*% alpha ) / (1 + exp(X %*% alpha ))^2
    w2_beta <- exp(X %*% alpha + beta) / (1 + exp(X %*% alpha + beta))^2

    # score for alpha and beta
    u_alpha <- t(X) %*% (Y - T0S * w1_alpha - T1S * w1_beta)
    u_beta <- sum(T1S * (Y - w1_beta))
    u <- rbind(u_alpha, u_beta)
    rownames(u) <- NULL

    # infomation matrix
    FI <- matrix(NA, d + 1, d + 1)
    D <- diag(c(T0S * w2_alpha + T1S * w2_beta))
    FI[1:d, 1:d] <- t(X) %*% D %*% X
    FI[d + 1, d + 1] <- sum(T1S * w2_beta)
    FI[1:d, d + 1] <- rowSums(t(X) %*% diag(c(T1S * w2_beta)))
    FI[d + 1, 1:d] <- rowSums(t(X) %*% diag(c(T1S * w2_beta)))
    theta <- tryCatch({
      c(alpha, beta) + c(solve(FI) %*% u)
    }, error = function(e){
      message("system is computationally singular")
      c(solve(FI +  0.01 * diag(d + 1)) %*%
          (FI %*% c(alpha, beta) + u))
    })

    if(max(abs(theta - theta_old)) < converged){
      break
    }else{
      theta_old <- theta
    }

  }

  return(list(alpha = theta[1:d], beta = theta[d + 1]))

}

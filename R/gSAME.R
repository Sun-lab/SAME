#' Association analysis using gene-level association analysis.
#'
#' @param Y The response variable. Could be continuous or binary.
#' @param X The design matrix. Intercept included.
#' @param O A matrix for the observed somatic mutation.
#' @param D A matrix for the total read-depth.
#' @param A A matrix for the number of alternative number matrix.
#' @param out_type The outcome type, "C" for continous, "D" for dichotomous. Default is "C".
#' @param theta_init The initail values of the parameters. Can be NULL.
#' @param mix_4bb A dataframe indicating the parameters of four beta-binomial distributions depending on the values of the observed somatic mutaton and the true somatic mutation when the read-depth is high.
#' @param null Logical. Indicating the estimation using EM algorim under the null hypothesis or not. The default is FALSE.
#' @param d0 The minimum of the total read-depth for obtaining the observed somatic mutation value. The default value is 20.
#' @param min_altcount The mimimum of the number of alternative reads that the somatic mutation could acutally orrur. The default value is 1.
#' @param spesen A dataframe specifying the specificity and sensitivity for all the somatic mutations.
#' @param mix_4bb A dataframe indicating the parameters of four beta-binomial distributions depending on the values of the observed somatic mutaton and the true somatic mutation when the read-depth is high.
#' @param mix_2bb A dataframe indicating the parameters of two beta-binomial distributions depending on the true value of the somatic mutation when the read-depth is low.
#' @param maxIT The maximal number of the EM iteration times. Default is 200.
#' @param converged The tolerance for the convergence. Default is 1e-6.
#' @param sig A matrix with the likelihood of O, D, A conditioning on the true value of the somatic mutation (S=0/1) for all the samples. Default is NULL.
#'
#' @return A list containing the output of the EM algorithm.
#' \item{Theta} A matrix including the estimators of the parameters in all the iterations.
#' \item{theta} The estimator of the parameters.
#' \item{LogLik} A vector of the log-likelihood values in all the iterations.
#' \item{logLik} The final log-likelihood.
#' \item{it} The number of the iteration times for convergence
#'
gSAME <- function(Y, X, O, D, A, out_type = "C",
                  theta_init, d0 = 20, null = FALSE,
                  sig = NULL, spesen = NULL, mix_2bb = NULL,
                  mix_4bb = NULL, min_altcount = 1,
                  maxIt = 200, converged = 1e-6, reEst = 1, traceIt = 0, ...){

  # X confounder matrix: n by d
  # O, D ,A matrices on gene level: n by p

  # compute likelihood of O, D, A if not given
  if(is.null(sig)){
    sig <- sigSAME(O, D, A, spesen = spesen, mix_4bb = mix_4bb,
                   mix_2bb = mix_2bb, min_altcount = min_altcount)
  }


  # estimators of confounder coefficients and noise level under null
  if(out_type == "C"){
    alpha_H0 <- solve( t(X) %*% X ) %*% t(X) %*% Y
    sigma2_H0 <- mean( (Y - X %*% alpha_H0)^2 )
  }else if(out_type == "D"){
    fit <- glm(Y ~ X - 1, family = "binomial")
    alpha_H0 <- coef(fit)
    sigma2_H0 <- 1
  }

  if(missing(theta_init)){
    theta_init <- c(alpha_H0, 0, sigma2_H0, 0.5)
  }

  theta_old <- theta_init
  theta <- theta_old

  n <- length(Y)
  d <- ncol(X)
  p <- ncol(O)

  Og <- (rowSums(O)>0) * 1

  Theta <- matrix(NA, maxIt, d + 3)

  LogLik = rep(NA, maxIt)

  for(it in 1:maxIt){

    if(traceIt){
      print(sprintf("%d: %s", it, date()))
    }

    alpha <- theta_old[1:d]
    beta <- ifelse(null,0,theta_old[d+1])
    sigma2 <- theta_old[d+2]
    rho0 <- theta_old[d+3]
    rho1 <- 1 - rho0

    ## E step
    # estimate class probability
    Mx <- matrix(NA, nrow = n, ncol = 2)

    for(i in 1:n){

      if(out_type == "C"){
        t01 <- rho0 * dnorm(Y[i], mean = sum(X[i,] * alpha),
                            sd = sqrt(sigma2))
        t11 <- rho1 * dnorm(Y[i],mean = sum(X[i,] * alpha) + beta,
                            sd = sqrt(sigma2))

      }else if(out_type == "D"){
        eta1 <- X %*% alpha
        eta2 <- X %*% alpha + beta
        liab1 <- exp(eta1)/(1 + exp(eta1) )
        liab2 <- exp(eta2)/(1 + exp(eta2) )
        t01 <- rho0 * dbinom(Y[i], 1, liab1[i])
        t11 <- rho1 * dbinom(Y[i], 1, liab2[i])
      }

      t02 <- sig[i, 1]
      t12 <- sig[i, 2]
      Mx[i,1] <- t01 * t02
      Mx[i,2] <- t11 * t12

    }

    T0S <- Mx[,1]/(Mx[,1]+Mx[,2])
    T1S <- Mx[,2]/(Mx[,1]+Mx[,2])

    # may cause unexpected result !!! Important!
    # T0S[is.na(T0S)] <- 1
    # T1S[is.na(T1S)] <- 0

    # for missing values; using the observed mutation as the true value
    # updated on 2018/07/23
    idna0 <- which(is.na(T0S))
    idna1 <- which(is.na(T1S))

    T1S[idna1] <- Og[idna1]
    T0S[idna0] <- 1 - T1S[idna0]


    Mx2 <- Mx[rowSums(Mx) != 0, ]

    LogLik[it]  <-  sum(log(Mx2[, 1] + Mx2[, 2]))


    if(reEst == 0){
      if(traceIt){
        print("reEst=0")
      }

      break
    }

    ## M step
    rho0 <- mean(T0S)
    rho0 <- min(rho0, 0.99)
    rho0 <- max(rho0, 0.01)

    rho1 <- 1 - rho0

    # M_tmp <- T1S %*% t(T1S) / sum(T1S)

    if(!null){
      if(out_type == "C"){
        M_tmp <- T1S %*% t(T1S) / sum(T1S)
        alpha <- solve(t(X) %*% (M_tmp - diag(n)) %*% X) %*% t(X) %*% (M_tmp - diag(n)) %*% Y
        alpha <- c(alpha)
        beta <- (sum(T1S * (Y - X %*% alpha)))/sum(T1S)
        r1 <- (Y - X %*% alpha)^2
        r2 <- (Y - X %*% alpha - beta)^2
        sigma2 <- sum(r1 * T0S) / n + sum(r2 * T1S) / n
      }else if(out_type == "D"){
        tmp <- fcmml(Y, X, alpha, beta, T0S, T1S)
        alpha <- tmp$alpha
        beta <- tmp$beta
        sigma2 <- 1
      }

    }else{
      alpha <- alpha_H0
      beta <- 0
      sigma2 <- sigma2_H0
    }


    theta <- c(alpha, beta, sigma2, rho0)

    Theta[it,] <- theta

    if(max(abs(theta - theta_old)) < converged){
      break
    }else{
      theta_old <- theta
    }

    if(traceIt){
      print(theta)
    }
  }

  names(theta) <-  c(paste("alpha",1:d,sep=""), "beta", "sigma2", "rho0")
  colnames(Theta) <- names(theta)
  mix1 <- list(Theta = Theta[1:it,], theta = theta,
               LogLik = LogLik[1:it], logLik = LogLik[it], it = it)
  mix1

}

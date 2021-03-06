#' Association analysis using mutation-level association analysis.
#'
#' @param Y The response variable. Could be continuous or binary.
#' @param X The design matrix. Intercept included.
#' @param O A vector for the observed somatic mutation.
#' @param D A vector for the total read-depth.
#' @param A A vector for the number of alternative number matrix.
#' @param out_type The outcome type, "C" for continous, "D" for dichotomous. Default is "C".
#' @param theta_init The initail values of the parameters. Can be NULL.
#' @param mix_4bb A dataframe indicating the parameters of four beta-binomial distributions depending on the values of the observed somatic mutaton and the true somatic mutation when the read-depth is high.
#' @param null Logical. Indicating the estimation using EM algorim under the null hypothesis or not. The default is FALSE.
#' @param d0 The minimum of the total read-depth for obtaining the observed somatic mutation value. The default value is 20.
#' @param gamm0 The specificity of the somatic mutation. Default is 1.
#' @param gamm1 The sensitivity of the somatic mutation. Default is 1.
#' @param bounds Some parameters for the bounds in the EM algorithm. Can be NULL.
#' @param maxIT The maximal number of the EM iteration times. Default is 200.
#' @param converged The tolerance for the convergence. Default is 1e-6.
#'
#' @return A list containing the output of the EM algorithm.
#' \item{Theta} A matrix including the estimators of the parameters in all the iterations.
#' \item{theta} The estimator of the parameters.
#' \item{LogLik} A vector of the log-likelihood values in all the iterations.
#' \item{logLik} The final log-likelihood.
#' \item{it} The number of the iteration times for convergence.


mSAME <- function(Y, X, O, D, A, out_type = "C", theta_init,
                 mix_4bb, null = FALSE, d0 = 20,
                 gamma0 = 1, gamma1 = 1, bounds = NULL, maxIt = 200,
                 converged = 1e-6, reEst = 1, traceIt = 0, ...){

  # bounds for parameters updated in the EM algorithm
  if(is.null(bounds)){
    bounds <- c(1e-4, 1e-4, 0.09, 0.01, 0.99)
  }

  min_pi <- bounds[1] # min pi for estimating bb parameters
  min_varphi <- bounds[2] # min varphi for estimating bb parameters
  pi0_upl <- bounds[3] # upper limit for pi0
  pi1_lwl <- bounds[4] # lower limit for pi1
  rho0_upl <- bounds[5] # upper limit for mixture proportion

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
    theta_init <- c(alpha_H0, 0, sigma2_H0, 0.5,
                    0.001, 0.005, 0.15, 0.40)
  }

  # the beta-binomal parameters of 4 mixtures (O = 0/1;S = 0/1)
  # default values are estimated from 3392 mutations
  if(missing(mix_4bb)){
    mix_4bb <- data.frame(type = c("TN", "FN", "FP", "TP"),
                          label = c(0, 2, 1, 3),
                          O = c(0, 0, 1, 1), S = c(0, 1, 0, 1),
                          pi = c(0.001029878, 0.036349845,
                                 0.046854682, 0.413395033),
                          varphi = c(0.0004398703, 0.3993310397,
                                     0.0001000000, 0.1497774366))
  }

  theta_old <- theta_init
  theta <- theta_old

  n <- length(Y)
  d <- ncol(X)
  Theta <- matrix(NA, maxIt, d + 7)

  LogLik <- rep(NA, maxIt)

  sam1 <- which(D < d0) # samIDs with low read-depth
  nTotal <- D[sam1]
  nA <- A[sam1]


  for(it in 1:maxIt){

    if(traceIt){
      print(sprintf("%d: %s", it, date()))
    }

    alpha <- theta_old[1:d]
    beta <- ifelse(null, 0, theta_old[d+1])
    sigma2 <- theta_old[d+2]
    rho0 <- theta_old[d+3]
    pi0 <- theta_old[d+4]
    varphi0 <- theta_old[d+5]
    pi1 <- theta_old[d+6]
    varphi1 <- theta_old[d+7]


    ## E step
    # estimate class probability
    Mx <- matrix(NA, nrow = n, ncol = 2)

    for(i in 1:n){

      if(out_type == "C"){
        t1 <- rho0 * dnorm(Y[i], mean = sum(X[i,] * alpha),
                           sd = sqrt(sigma2))
        t3 <- (1 - rho0) * dnorm(Y[i],mean = sum(X[i,] * alpha) + beta,
                               sd = sqrt(sigma2))

      }else if(out_type == "D"){
        eta1 <- X %*% alpha
        eta2 <- X %*% alpha + beta
        liab1 <- exp(eta1)/(1 + exp(eta1) )
        liab2 <- exp(eta2)/(1 + exp(eta2) )
        t1 <- rho0 * dbinom(Y[i], 1, liab1[i])
        t3 <- (1 - rho0) * dbinom(Y[i], 1, liab2[i])
      }

      t2 <- ifelse(i %in% sam1, dbetabinom(A[i], size = D[i], prob = pi0, rho = varphi0),
                   (O[i]==0) * gamma0 * dbetabinom(A[i], size = D[i],
                                                   prob = mix_4bb$pi[1],
                                                   rho = mix_4bb$varphi[1]) +
                   (O[i]==1) * (1 - gamma0) * dbetabinom(A[i], size = D[i],
                                                         prob = mix_4bb$pi[2],
                                                         rho = mix_4bb$varphi[2]))

      t4 <- ifelse(i %in% sam1, dbetabinom(A[i], size = D[i], prob = pi1, rho = varphi1),
                   (O[i]==1) * gamma1 * dbetabinom(A[i], size = D[i],
                                                   prob = mix_4bb$pi[4],
                                                   rho = mix_4bb$varphi[4]) +
                   (O[i]==0) * (1 - gamma1) * dbetabinom(A[i], size = D[i],
                                                         prob = mix_4bb$pi[3],
                                                         rho = mix_4bb$varphi[3]))

      Mx[i,1] <- t1 * t2
      Mx[i,2] <- t3 * t4

    }

    T0S <- Mx[,1] / (Mx[,1]+Mx[,2])
    T1S <- Mx[,2] / (Mx[,1]+Mx[,2])

    # may cause unexpected result !!! Important!
    # T0S[is.na(T0S)] <- 1
    # T1S[is.na(T1S)] <- 0

    # for missing values; using the observed mutation as the true value
    # updated on 2018/07/23
    idna0 <- which(is.na(T0S))
    idna1 <- which(is.na(T1S))

    T1S[idna1] <- O[idna1]
    T0S[idna0] <- 1 - T1S[idna0]


    Mx2 <- Mx[rowSums(Mx) != 0, ]
    LogLik[it]  <-  sum(log(Mx2[,1] + Mx2[,2]))

    if(reEst == 0){
      if(traceIt){
        print("reEst=0")
      }

      break
    }

    ## M step
    rho0 <- mean(T0S)
    rho0 <- min(rho0, rho0_upl)
    rho0 <- max(rho0, 0.01)

    log0 <- capture.output({
      fit0 <- mle.bb(nA, nTotal, pi0, varphi0,
                     min_pi, min_varphi, ws = T0S[sam1])
    })

    pi0 <- min(fit0$pi, pi0_upl)
    varphi0 <- fit0$rho

    log1 <- capture.output({
      fit1 <- mle.bb(nA, nTotal, pi1, varphi1,
                     min_pi, min_varphi, ws = T1S[sam1])
    })

    pi1 <- max(fit1$pi, pi1_lwl)
    varphi1 <- fit1$rho

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

    theta <- c(alpha, beta, sigma2, rho0, pi0, varphi0, pi1, varphi1)

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

  names(theta) <-  c(paste("alpha",1:d,sep=""),"beta","sigma2","rho0",
                     "pi0", "varphi0", "pi1", "varphi1")
  colnames(Theta) <- names(theta)
  mix1 <- list(Theta = Theta[1:it, ], theta = theta, LogLik = LogLik[1:it], logLik = LogLik[it], it = it)
  mix1

}

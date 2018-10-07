
#' Estimate the specificity and sensitivity of a somatic mutation based on its read-count data.
#'
#' @param A A vector for the number of alternative reads for the somatic mutation.
#' @param D A vector for the read-depths for the somatic mutation.
#' @param O A vector for the observed somatic mutation values.
#' @param pi The mean parameters for the beta-binomial distributions when the true somatic mutation value equals 0 or 1.
#' @param varphi The over-dispersion parameters for the beta-binomial distributions when the true somatic mutation value equals 0 or 1.
#'
#' @return The estimated specificity and sensitivity for the somatic mutaiton.




estss <- function(A, D, O, pi, varphi){

  pi0 <- pi[1]
  varphi0 <- varphi[1]
  pi1 <- pi[2]
  varphi1 <- varphi[2]

  rho0 <- 1 - mean(O)
  rho1 <- 1 - rho0

  nA <- A
  nTotal <- D
  # estimate class probability
  logL_Mx <- matrix(NA, nrow = length(nTotal), ncol = 2)

  for(i in 1:length(nTotal)){
    logL_Mx[i,1] <- dbetabinom(nA[i], size = nTotal[i], prob = pi0,
                              rho = varphi0, log = TRUE)
    logL_Mx[i,2] <- dbetabinom(nA[i], size = nTotal[i], prob = pi1,
                              rho = varphi1, log = TRUE)
  }

  logL_Mx[,1] <- logL_Mx[,1] + log(rho0)
  logL_Mx[,2] <- logL_Mx[,2] + log(rho1)

  # calculate posterior prob while keeping numbers in log scale
  postPs <- 1 / (1 + exp(logL_Mx[,2] - logL_Mx[,1]))
  pp <- cbind(postPs, 1 - postPs)

  O_true <- (postPs <= 0.5) * 1

  tp <- sum(O_true == 1 & O == 1)
  tn <- sum(O_true == 0 & O == 0)

  fp <- sum(O_true == 0 & O == 1)
  fn <- sum(O_true == 1 & O == 0)

  if(tp + fn > 0){
    gamma1 <- tp / (tp + fn)
  }else if(tp + fn == 0){
    gamma1 <- 1
  }
  gamma0 <- tn / (tn + fp)
  est <- c(gamma0, gamma1)
  names(est) <- c("spe", "sen")
  est

}

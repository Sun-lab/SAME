#' Compute the likelihood for the observed somatic mutation and the read-depth data for gSAME.
#' It can significantly reduce computational time for assication analysis
#' between one gene-level mutation and lots of outcomes.
#'
#' @param O A matrix for the observed somatic mutation.
#' @param D A matrix for the total read-depth.
#' @param A A matrix for the number of alternative number matrix.
#' @param d0 The minimum of the total read-depth for obtaining the observed somatic mutation value. The default value is 20.
#' @param min_altcount The mimimum of the number of alternative reads that the somatic mutation could acutally orrur. The default value is 1.
#' @param spesen A dataframe specifying the specificity and sensitivity for all the somatic mutations.
#' @param mix_4bb A dataframe indicating the parameters of four beta-binomial distributions depending on the values of the observed somatic mutaton and the true somatic mutation when the read-depth is high.
#' @param mix_2bb A dataframe indicating the parameters of two beta-binomial distributions depending on the true value of the somatic mutation when the read-depth is low.
#'
#' @return A matrix with the likelihood of O, D, A conditioning on the true value of the somatic mutation (S=0/1) for all the samples.









sigSAME <- function(O, D, A, d0 = 20, min_altcount = 1,
                    spesen = NULL, mix_4bb = NULL, mix_2bb = NULL){
  n <- nrow(O)
  p <- ncol(O)

  # specificity and sensitivity
  if(is.null(spesen)){
    spesen <- data.frame(spe = rep(1,p),sen = rep(1,p))
  }
  # which mutations have calling errors
  index2 <- which(rowSums(spesen) < 2)


  # parameters for 4 bb mixtures for high read-depth
  if(is.null(mix_4bb)){
    mix_4bb <- data.frame(type = c("TN","FN","FP","TP"),
                          label = c(0,2,1,3),
                          O = c(0,0,1,1), S = c(0,1,0,1),
                          pi = c(0.0009500459, 0.0019825522, 0.1178939318, 0.3206793574),
                          varphi = c(0.0006210529, 0.3456750569, 0.0001000000, 0.1018129448))
  }

  # parameters for 2 bb mixtures for low read-depth
  nA <- A[D < d0]
  nTotal <- D[D < d0]
  if(is.null(mix_2bb)){
    if(length(nA) > 0){
      log01 <- capture.output({
        fit01 <- bb.mix(nA,nTotal)
      })
      mix_2bb <- fit01$theta[1:4]
    }else{
      mix_2bb <- c(0.001,0.005,0.15,0.40)
    }
  }
  pi0 <- mix_2bb[1]
  varphi0 <- mix_2bb[2]
  pi1 <- mix_2bb[3]
  varphi1 <- mix_2bb[4]

  df <- data.frame(O = O, D = D, A = A)
  # df <- cbind(O, D, A)

  freq <- unname(colMeans(O))

  saveInfo <- function(x){
    # x <- unlist(x)

    Oi <- x[1:p]
    Di <- x[(1 + p):(2 * p)]
    Ai <- x[(1 + 2 * p):(3 * p)]

    # function to compute the likelihood on the mutation level
    fmut <- function(A, D, O, S, spe, sen){
      if(D < d0){
        if(S == 0){
          res <- dbetabinom(A, size = D, prob = pi0, rho = varphi0)
        }else if(S == 1){
          res <- dbetabinom(A, size = D, prob = pi1, rho = varphi1)
        }
      }else{
        prob <- ifelse(S == 0, 1 - spe, sen)
        index_tmp <- match(O + 2 * S, mix_4bb$label)
        res <- ifelse(O == 0,1 - prob, prob) *
          dbetabinom(A, size = D, prob = mix_4bb$pi[index_tmp],
                     rho = mix_4bb$varphi[index_tmp])
      }
      res
    }




    # case 1. gene-mut = 0 =====
    prob_mut0 <- Vectorize(fmut)(A = Ai, D = Di, O = Oi, S = 0,
                                 spe = spesen[, 1], sen = spesen[, 2])
    t02 <- prod(prob_mut0)


    # case 2. gene-mut = 1 ======

    # compute the mut-level likelihood of A,D,O conditional on S
    # assume the mutation can orrur only when A_{ij} >= 1,

    # further reduce the computation burden by only considering
    # the mutations with calling errors (spe != 1 or sen ! = 1)

    # which mutations may have true value 1
    index1 <- which(Ai >= min_altcount) # length k

    if(length(index1) != 0){

      # which mutations may have true value 1 and have calling errors
      index31 <- index1[index1 %in% index2]

      # which mutations may have true value 1 and NOT have calling errors
      index32 <- index1[!index1 %in% index2]

      if(length(index31) != 0){
        grid <- as.matrix(expand.grid(rep(list(0:1),length(index31))))
        setS <- matrix(0, 2^(length(index31)), ncol(O))
        setS[, index31] <- grid
        if(length(index32) >= 1){
          setS[, index32] <-  matrix(rep(Oi[index32], 2^(length(index31))),
                                     ncol = 2^(length(index31)), byrow = T)
        }

        rs <- rowSums(setS)
         setS <- setS[rs != 0, , drop = FALSE]

        prob1 <- apply(setS, 1, function(x) abs(prod((x==0) - freq)))

        prob1 <- prob1/sum(prob1)

        # p by 2^k-1
        prob2 <- apply(setS,1,function(x){
          Vectorize(fmut)(A = Ai, D = Di, O = Oi, S = x,
                          spe = spesen[, 1], sen = spesen[, 2])
        })

        prob2 <- as.matrix(prob2)
        v12 <- prob1 * apply(prob2, 2, prod)

        t12 <- sum(v12)

      }else{
        # mutations either Ai == 0 or no calling error
        Si <- rep(0, p)
        Si[index32] <- Oi[index32]
        if(sum(Si) != 0 ){ # the true mut != 0
          prob_mut1 <- Vectorize(fmut)(A = Ai, D = Di, O = Oi, S = Si,
                                       spe = spesen[, 1], sen = spesen[, 2])
          t12 <- prod(prob_mut1)
        }else{ # the true mut == 0
          t12 <- 0
        }

      }
    }else{
      t12 <- 0
    }

    return(c(t02 = t02, t12 = t12))
  }
  sig <- t(apply(df, 1, saveInfo))
  sig
}

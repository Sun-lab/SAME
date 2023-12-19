
library(devtools)
library(VGAM)
library(MASS)
library(basic) # install_github("Sun-lab/basic/basic")
library(SAME) # install_github("Sun-lab/SAME")


# mSAME -----

# import the sample mutation, total read-depth, alternative read-depth data
mutM <- read.table("mut_dat.txt", sep = "\t", header = TRUE)
Depth_tumor <- read.table("depth_dat.txt", sep = "\t", header = TRUE)
Alter_tumor <- read.table("alter_dat.txt", sep = "\t", header = TRUE)

# import the response and adjusting covariates data
dat_reg <- read.table("dat_reg.txt", sep = "\t", header = TRUE)
Y <- dat_reg$Y
X <- cbind(1, dat_reg$X)

# estimate the parameters of the mixed beta binomial distributions 
d0 <- 20

D0 <- Depth_tumor[mutM==0 & Depth_tumor >= d0] 
A0 <- Alter_tumor[mutM==0 & Depth_tumor >= d0]
D1 <- Depth_tumor[mutM==1 & Depth_tumor >= d0] 
A1 <- Alter_tumor[mutM==1 & Depth_tumor >= d0]

fit0_mut <- bb.mix2(A0,D0) 
fit1_mut <- bb.mix2(A1,D1)

mix_4bb_hat <- data.frame(type = c("TN","FN","FP","TP"),
                          label = c(0,2,1,3),
                          O = c(0,0,1,1),
                          S = c(0,1,0,1),
                          pi = NA,
                          varphi = NA)

mix_4bb_hat$pi <- c(fit0_mut$theta[1],fit0_mut$theta[3],
                    fit1_mut$theta[1],fit1_mut$theta[3])

mix_4bb_hat$varphi <- c(fit0_mut$theta[2],fit0_mut$theta[4],
                        fit1_mut$theta[2],fit1_mut$theta[4])

# estimate specificity & sensitivity
fit_ss0 <- mle.bb(A0, D0)
fit_ss1 <- mle.bb(A1, D1)
estss_pi <- c(fit_ss0$pi, fit_ss1$pi)
estss_varphi <- c(fit_ss0$rho, fit_ss1$rho)

p <- ncol(mutM)

spesen_est <- matrix(NA, p, 2)
for(j in 1:p){
  O <- mutM[, j] 
  D <- Depth_tumor[, j]
  A <- Alter_tumor[, j]
  spesen_est[j, ] <- SAME::estss(A, D, O, estss_pi, estss_varphi)
}


# compute p-value for each of the 20 mutations in mSAME model
pval <- c()
for(j in 1:p){
  
  # fit the model
  fit <- mSAME(Y, X, O, D, A, gamma0 = spesen_est[j, 1], 
               gamma1 = spesen_est[j, 2],
               mix_4bb = mix_4bb_hat)

  # fit the model under H0
  fit0 <- mSAME(Y, X, O, D, A, null = TRUE, 
                gamma0 = spesen_est[j, 1], 
                gamma1 = spesen_est[j, 2], 
                mix_4bb = mix_4bb_hat)
  t1 <- 2 * (fit$logLik - fit0$logLik) # LRT statistic 
  pval[j] <- exp(pchisq(t1, df = 1, log.p = T, lower.tail = F))
  
}

# gSAME ----
# use the same data above 
# but consider one gene-level mutation with 20 individual mutations


O <- mutM
D <- Depth_tumor
A <- Alter_tumor

# estimate the extra parameters of the mixed beta binomial in gSAME
nA <- A[D < d0]
nTotal <- D[D < d0]
if(length(nA) > 0){
  log01 <- capture.output({
    fit01 <- bb.mix2(nA, nTotal)
  })
  mix_2bb <- fit01$theta[1:4]
}else{
  mix_2bb <- c(0.001, 0.005, 0.15, 0.40)
}

# Compute likelihood for the observed somatic mutation and the read-depth data
# to reduce the computational cost in gSAME
sig <- sigSAME(O, D, A, spesen = spesen_est, mix_4bb = mix_4bb_hat)

fit_gSAME <- gSAME(Y, X, O, D, A, spesen = spesen_est, 
                  mix_4bb = mix_4bb_hat,
                  sig = sig, mix_2bb = mix_2bb)

fit0_gSAME <- gSAME(Y, X, O, D, A, null = TRUE, spesen = spesen_est, 
                    mix_4bb = mix_4bb_hat,
                    sig = sig, mix_2bb = mix_2bb)

t2 <- 2*(fit_gSAME$logLik - fit0_gSAME$logLik)
pval_gSAME <- exp(pchisq(t2, df = 1, log.p = T, lower.tail = F))














library("ashr")
library("mashr")
set.seed(1)
heri = "0.5"#heritability of the simulated data
causal = "1000" #the number of causal SNPs in the simulated data
mash_BETA <- read.csv(paste0("beta/beta.heri.",heri, ".causal.",causal, ".ampli.txt", header = TRUE, sep = " ")
mash_SE <- read.csv(paste0("se/se.heri.",heri,".causal.",causal".ampli.txt"), header = TRUE, sep = " ")
mash_BETA <- as.matrix(mash_BETA) ; mash_SE <- as.matrix(mash_SE)
data <- mash_set_data(mash_BETA, mash_SE, zero_Shat_reset = .Machine$double.eps)
U.c = cov_canonical(data)
n <- 30
sigma <- matrix(numeric(n^2), n, n)
for (i in 1:n) {
  for (j in 1:n) {
    if(i==j){
      sigma[i,j] <- 1# the variance is set as equal: tbd MARK
    }else{
      sigma[i,j] <- 0.9^(abs(i-j))#simulate the decreasing covariance: possible functions: 0.9^{age_i - age_j}; 
    }
  }
}
tmp <- data.frame("mat","weight")
  # fit mash model 
  U.pca <- cov_pca(data, 5)
  U.ed <- cov_ed(data, U.pca) #subset = strong)

  m = mash(data, Ulist= c(U.c, U.ed)
  pos_mean <- get_pm(m)
  pos_sd <- get_psd(m)
  lfsr <- get_lfsr(m)
  weight <- get_estimated_pi(m)
  Sys.time()
  write.table(pos_mean, "pos_estimates/pos_m.0.5.heri.2000.ampli.txt")
  write.table(pos_sd, "pos_estimates/pos_sd.0.5.heri.2000.ampli.txt")
  write.table(lfsr, "lfsr/lfsr.0.7.heri.2000.ampli.txt")
  write.table(weight,"weight/weight.0.5.heri.2000.ampli.txt")
  png(filename = "heatmap.U.ed.0.5.2000.png")
  heatmap(U.ed)
  dev.off()
  barplot(weight,las = 2)
  mash_plot_meta(m,get_significant_results(m)[1])

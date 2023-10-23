# Load the required libraries
library("ashr")
library("mashr")

# Set a seed for reproducibility
set.seed(1)

# Define the heritability and the number of causal SNPs for the simulated data
heri = "0.5"
causal = "1000"

# Read in the BETA and SE data from specified files
mash_BETA <- read.csv(paste0("beta/beta.heri.",heri, ".causal.",causal, ".ampli.txt"), header = TRUE, sep = " ")
mash_SE <- read.csv(paste0("se/se.heri.",heri,".causal.",causal,".ampli.txt"), header = TRUE, sep = " ")

# Convert the data frames to matrices
mash_BETA <- as.matrix(mash_BETA)
mash_SE <- as.matrix(mash_SE)

# Set up the data for mash analysis
data <- mash_set_data(mash_BETA, mash_SE, zero_Shat_reset = .Machine$double.eps)

# Get canonical covariance matrices
U.c = cov_canonical(data)

# Define the dimension and initialize the sigma matrix
n <- 30
sigma <- matrix(numeric(n^2), n, n)

# Fill in the sigma matrix
for (i in 1:n) {
  for (j in 1:n) {
    if(i==j){
      sigma[i,j] <- 1  # Set variance on the diagonal
    }else{
      sigma[i,j] <- 0.9^(abs(i-j))  # Set covariance for off-diagonal elements based on distance
    }
  }
}

# Initialize a data frame to store matrices and weights
tmp <- data.frame("mat","weight")

# Fit the mash model
U.pca <- cov_pca(data, 5)  # Principal component analysis on the data
U.ed <- cov_ed(data, U.pca)  # Empirical Bayes decomposition

m = mash(data, Ulist= c(U.c, U.ed))  # Fit the mash model with the specified covariance matrices
pos_mean <- get_pm(m)  # Get posterior mean
pos_sd <- get_psd(m)  # Get posterior standard deviation
lfsr <- get_lfsr(m)  # Get local false sign rate
weight <- get_estimated_pi(m)  # Get estimated mixture weights

# Write the results to files
write.table(pos_mean, "pos_estimates/pos_m.0.5.heri.2000.ampli.txt")
write.table(pos_sd, "pos_estimates/pos_sd.0.5.heri.2000.ampli.txt")
write.table(lfsr, "lfsr/lfsr.0.7.heri.2000.ampli.txt")
write.table(weight,"weight/weight.0.5.heri.2000.ampli.txt")

# Create a heatmap of the empirical Bayes decomposition matrix
png(filename = "heatmap.U.ed.0.5.2000.png")
heatmap(U.ed)
dev.off()

# Create a barplot of the estimated mixture weights
barplot(weight,las = 2)

# Plot meta-analysis results for significant findings
mash_plot_meta(m,get_sign

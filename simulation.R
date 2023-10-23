# This script is utilized for generating simulated phenotype data.
library("MASS")  # Load the MASS library for generating multivariate normal random variates

order <- as.integer(args[1])  # Convert the first argument to integer
snp_num <- 2000  # Define the number of causal SNPs
set.seed(order)  # Set the random seed for reproducibility

# Simulating age data using a multinomial distribution
num_samples <- 300000  # Define the number of samples
age_probs <- rep(1/30, 30)  # Assuming 30 possible age values
simulated_ages <- sample(1:30, size = num_samples, replace = TRUE, prob = age_probs)  # Simulate ages
age_indices <- list()
for (age_value in 1:30) {
  indices <- which(simulated_ages == age_value)
  age_indices[[paste("id", age_value, sep = "_")]] <- indices  # Store indices for each age group
}

# Reading and processing genotype matrix
snp_freqs <- read.csv("maf_sample_20k.txt", colClasses = 'numeric')  # Read genotype data
snp_freqs <- snp_freqs$x
load("simulation_matrix_k_5k.RData")  # Load genotype matrix data from provided file
snp_freqs <- snp_freqs[1:snp_num]  # Select causal SNPs
genotype_matrix_i <- genotype_matrix_k[,(1:snp_num)]
dim(genotype_matrix_i)  # Display dimensions of genotype matrix

# Generating Beta (effect size) matrix
print("create matrix of Betas sampled from prespecified covariance matrices")
### PRE-SPECIFIED COV MATRICES ###
n <- 30
mu <- rep(0,n)  # Mean of the effect sizes is set to 0
# Create covariance matrix
sigma <- matrix(numeric(n^2), n, n)
for (i in 1:n) {
  for (j in 1:n) {
    if(i==j){
      sigma[i,j] <- 1.1^i  # Set variance (diagonal of covariance matrix)
    }else{
      sigma[i,j] <- 0.9^(abs(i-j))  # Set covariance for off-diagonal elements
    }
  }
}
Beta <- mvrnorm(snp_num, mu=mu, Sigma=sigma)  # Generate effect sizes using multivariate normal distribution

# Environmental effect computation
print("get vector of environmental effect")
E <- get_environment(0.5)  # Get environmental effect vector with heritability set to 0.5
length(E)  # Display length of environmental effect vector

# Cleanup
rm(snp_freqs)

# Generating phenotype vector
print("get phenotype vector")
genotype_matrix_i <- t(genotype_matrix_i)  # Transpose genotype matrix
for (age_value in 1:30) {
  genotype_matrix_i[,age_indices[[paste("id", age_value, sep = "_")]]] <- genotype_matrix_i[,age_indices[[paste("id", age_value, sep = "_")]]] * Beta[,age_value]  # Multiply genotypes by effect sizes
}
length(genotype_matrix_i)
pheno <- colSums(genotype_matrix_i)+ E  # Sum across columns and add environmental effect
save(pheno,age_indices, file = "pheno.2000.0.5.RData")  # Save phenotype data and age indices to file
print(length(pheno))
rm(genotype_matrix_i)  # Cleanup
# Define the GWAS function to perform a linear regression analysis on genotype and phenotype data for a specified age group
GWAS <- function(genotype, pheno, E, age_value) {
  x <- genotype[age_indices[[paste("id", age_value, sep = "_")]]]  # Extract genotype data for the specified age group
  y <- pheno[age_indices[[paste("id", age_value, sep = "_")]]]  # Extract phenotype data for the specified age group
  e <- E[age_indices[[paste("id", age_value, sep = "_")]]]  # Extract environmental effect data for the specified age group
  model <- lm(y ~ x, na.action = na.omit)  # Fit a linear model while omitting missing values
  return(summary(model))  # Return the summary of the linear model
}
print("perform GWAS")

# Initialize lists and data frames to store beta values and standard errors
beta_value <- list()
se_value <- list()
mash_BETA <- data.frame(matrix(NA, nrow = num_samples, ncol = 30))  # Data frame to store beta values
mash_SE <- data.frame(matrix(NA, nrow = num_samples, ncol = 30))  # Data frame to store standard errors

# Loop through each age group and perform GWAS
for (age in 1:30) {
    print(paste("age is ", age))
    beta_value[[age]] <- numeric(5000)  # Initialize numeric vector to store beta values
    se_value[[age]] <- numeric(5000)  # Initialize numeric vector to store standard errors
    for (i in 1:length(colnames(genotype_matrix_k))) {
        gwas <- GWAS(genotype_matrix_k[,i], pheno, E, age)  # Perform GWAS for each SNP within the age group
        beta_value[[age]][i] <- gwas$coefficients[2, 1]  # Store beta value
        se_value[[age]][i] <- gwas$coefficients[2, 2]  # Store standard error
    }
    mash_BETA[,age] <- unlist(beta_value[[age]])  # Populate data frame with beta values
    mash_SE[,age] <- unlist(se_value[[age]])  # Populate data frame with standard errors
}

# Clean up by removing phenotype_matrix_k
rm(phenotype_matrix_k)

# Write beta values and standard errors to files
write.table(mash_BETA ,paste0("beta/beta.",".heri.0.5.causal.2000.ampli.txt") )  # Write beta values to file
write.table(mash_SE, paste0("se/se.",".heri.0.5.causal.2000.ampli.txt"))  # Write standard errors to file

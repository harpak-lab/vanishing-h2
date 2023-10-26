library("ashr")
library("mashr")
args = commandArgs(trailingOnly=TRUE)
trait <- args[1]
threshold <- "1e-5" 
## read in BETA and SE 
read.sum <- function(age){
  aar <- read.table(paste(trait, "/aar.", age, ".", trait, ".",threshold ,".glm.linear", sep = ""), header = TRUE, col.names = c("VAR", "BETA", "SE"))
  return(aar)
}

ages <- 40:69
datasets <- lapply(ages, function(age) read.sum(age))
names(datasets) <- paste("aar", ages, sep = "")
BETA_list <- lapply(datasets, function(dataset) dataset$BETA)
SE_list <- lapply(datasets, function(dataset) dataset$SE)
if(is.null(names(BETA_list)) || is.null(names(SE_list))) {
  stop("BETA_list and SE_list must be named lists.")
}
names(BETA_list)
BETA <- do.call(cbind, BETA_list)
SE <- do.call(cbind, SE_list)
head(BETA)
colnames(BETA) <- names(BETA_list)
colnames(SE) <- names(SE_list)
if("VAR" %in% colnames(BETA_list[[1]])) {
  rownames(BETA) <- BETA_list[[1]]$VAR
  rownames(SE) <- SE_list[[1]]$VAR
}


print("This is the original BETA and SE")
head(BETA)
dim(BETA)
head(SE)
dim(SE)

#Remove rows with NA 
na_rows <- rowSums(is.na(BETA) | is.na(SE)) > 0
BETA <- BETA[!na_rows, ]
SE <- SE[!na_rows, ]
print("This is the BETA and SE after filtering NAs")
dim(BETA)
dim(SE)
head(BETA)
head(SE)
header <- as.character(40:69)
pm_all <- data.frame(matrix(ncol = 30, nrow = 0)) ; colnames(pm_all) <- header
psd_all <- data.frame(matrix(ncol = 30, nrow = 0)) ; colnames(psd_all) <- header
lfsr_all <- data.frame(matrix(ncol = 30, nrow = 0)) ; colnames(lfsr_all) <- header
##mash
data = mash_set_data(BETA, SE)
U.pca <- cov_pca(data, 5)
U.ed <- cov_ed(data, U.pca)#data-driven matrices
U.c = cov_canonical(data)
# fit mash model
m = mash(data, Ulist= c(U.c, U.ed))
# mixture model
pos_mean <- get_pm(m)
pos_sd <- get_psd(m)
lfsr <- get_lfsr(m)

write.table(pos_mean, file = paste(trait,"/pos.mean.", threshold, ".",trait,".txt", sep = ""))
write.table(pos_sd, file = paste(trait, "/pos.sd.",threshold,".", trait, ".txt", sep = ""))
write.table(lfsr, file = paste(trait, "/lfsr.", threshold,".",trait, ".txt", sep = ""))

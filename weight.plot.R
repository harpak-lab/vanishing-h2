#The script is to plot the matrix weight for the simulated data(10 different datasets)
library(ggplot2)

heri <- "0.5"
causal <- "1000"

# Read files and add matrix column
df_names <- paste0("weight", 11:20)

for(i in seq_along(df_names)) {
  file_name <- paste0("weight/weight.", heri, ".heri.", causal, ".causal.", i, ".txt")
  df <- read.csv(file_name, header = TRUE, sep = " ", col.names = "weight")
  df$matrix <- row.names(df)
  assign(df_names[i], df)
}
# Combine weights from all dataframes into one matrix
tmp <- do.call(cbind, sapply(df_names, function(x) get(x)$weight, simplify = FALSE))
colnames(tmp) <- df_names
tmp <- as.data.frame(tmp)
# Compute average and standard error
tmp$se <- apply(tmp, 1, function(x) sd(x) / sqrt(length(x)))
tmp$mean <- rowMeans(tmp)
tmp$matrix <- row.names(weight1)
tmp
p <- ggplot(data = tmp, aes(x = matrix, y = mean)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  labs(
    title = "Matrix Weight simulation",
    subtitle = paste("heritability =", heri, "; causal SNPs =", causal),
    x = "matrix weight",
    y = "covariance matrix"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(filename = paste0("weight/matrix.weight.heri.", heri, ".causal.", causal, ".pdf"), plot = p)

# Set working directory
setwd("/home/ibg-4/Desktop/Rhome/moca_blue/mo_nom")
# Load required libraries
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("TFBSTools")
#library(TFBSTools)
#library(motifStack)
#library(universalmotif)
#library(JASPAR2020)
# Update the file name and path accordingly
jaspar_file <- "./out/rdf5_2At-Sb-Sl-ZmM0+S0_cwm-motifs.jaspar"
pfm <- read_jaspar(jaspar_file)
pwm_uni0<-convert_motifs(pfm
  , class = "TFBSTools-PWMatrix")
jaspar_motifs_plants <- getMatrixSet(JASPAR2020,
                                     opts = list(tax_group='plants',
                                                 matrixtype='PWM'))
###################################################
similarities_list <- vector("list", length = length(pwm_uni0))
####################################################
for (i in seq_along(pwm_uni0)) {
  similarities <- numeric(length(jaspar_motifs_plants))
  
  for (j in seq_along(jaspar_motifs_plants)) {
    similarities[j] <- PWMSimilarity(pwm_uni0[[i]]@profileMatrix,
                                     jaspar_motifs_plants[[j]]@profileMatrix,
                                     method = 'Pearson')
  }
    similarities_list[[i]] <- similarities
}
####################################################
max_values <- numeric(length(similarities_list))
max_positions <- numeric(length(similarities_list))
for (i in seq_along(similarities_list)) {
  max_values[i] <- max(similarities_list[[i]])
  max_positions[i] <- which.max(similarities_list[[i]])
}
max_table <- data.frame(Idx= c(1:length(pwm_uni0)),
                        MaxValue = max_values,
                        MaxPosition = max_positions)
#####################################################
pwm_list <- lapply(pwm_uni0, function(pwm) {
  pwm@name
})
pwm_df <- as.data.frame(do.call(rbind, pwm_list))
pwm_df$Idx <- c(1:nrow(pwm_df))
names(pwm_df)[1]<- "epm"
df0 <- merge(pwm_df, max_table, by= "Idx")
###
jaspar_motifs_plants_list<- as.data.frame(names(jaspar_motifs_plants))
jaspar_motifs_plants_list$MaxPosition <- c(1:nrow(jaspar_motifs_plants_list))
df1 <- merge(df0, jaspar_motifs_plants_list, by = "MaxPosition")
##
df1 <- df1[, !(names(df1) %in% c("MaxPosition", "Idx"))]
colnames(df1) <- c( "epm", "maxPearson", "TFBS_JASPAR2020")
#####################################################

#####################################################
random_values <- runif(100000)
set.seed(1749)
sampled_values <- sample(random_values, length((similarities_list[[1]])))

p_vals <- list()
for (i in seq_along(similarities_list)) {
#  set.seed(1749)
   sampled_values <- sample(random_values, length((similarities_list[[1]])))
   mean_value <- mean(sampled_values)
#  sd_value <- sd(sampled_values)
#  lower_bound <- (max(similarities_list[[i]])) - max(similarities_list[[i]]) * 0.05
#  upper_bound <- (max(similarities_list[[i]])) + max(similarities_list[[i]]) * 0.05
#  p_val <- pnorm(upper_bound, mean = mean_value, sd = sd_value) - pnorm(lower_bound, mean = mean_value, sd = sd_value)
  mean_value <- mean(sampled_values)
  sd_value <- sd(sampled_values)
  p_val <- pnorm(max(similarities_list[[i]]), mean = mean_value, sd = sd_value)
  p_vals[[i]] <- p_val
}

print(p_vals)

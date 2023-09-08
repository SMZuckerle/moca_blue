file_path <- "tissue_predictions_results.csv"
data <- read.csv(file_path)
#################################
########################################################################################
#set.seed(42)  # For reproducibility
N <- 1000  # Number of iterations
# Number of samples in the test sets
test_set_size <- 100
# Bootstrapping function
bootstrap <- function(data) {
  n <- nrow(data)
  sample_indices <- sample(1:n, size = test_set_size, replace = TRUE)
  return(data[sample_indices, ])
}
###########################################################################################
data$root_leaf_specificity <- ifelse(
  data$true_target_root == 0 & data$true_target_leaf == 1,
  "root_low_leaf_high",
  ifelse(
    data$true_target_root == 0 & data$true_target_leaf == 0,
    NA,
    "leaf_low_root_high"
  )
)
colnames(data)
data$root_leaf_correct <- ifelse(
  data$pred_root ==  data$true_target_leaf,
  "True", "False")
data$leaf_root_correct <- ifelse(
  data$pred_leaf ==  data$true_target_root,
  "True", "False")
split_data <- split(data, data$specie)
lapply(split_data, head)
##########################################################################################
arabidopsis_df <- split_data$arabidopsis
sbicolor_df <- split_data$sbicolor
solanum_df <- split_data$solanum
zea_df <- split_data$zea
###########################################################################################
specific_genes <- zea_df[zea_df$specific == "specific", ]
non_specific_genes <- zea_df[zea_df$specific == "non-specific", ]
specific_genes_root_0_leaf_1 <- specific_genes[specific_genes$root_leaf_specificity == "root_low_leaf_high", ]
specific_genes_root_1_leaf_0 <- specific_genes[specific_genes$root_leaf_specificity == "leaf_low_root_high", ]
###########################################################################################
####################################
results <- data.frame(
  iteration = integer(N),
  Accuracy_non_specific_leaf_leaf = double(N),
  Bccuracy_non_specific_leaf_root = double(N),
  Cccuracy_non_specific_root_root = double(N),
  Dccuracy_non_specific_root_leaf = double(N),
  Eccuracy_specific_leaf_leaf = double(N),
  Fccuracy_specific_leaf_root = double(N),
  Gccuracy_specific_root_root = double(N),
  Hccuracy_specific_root_leaf = double(N),
  Iccuracy_specific_genes_root_0_leaf_1_leaf_leaf = double(N),
  Jccuracy_specific_genes_root_0_leaf_1_leaf_root = double(N),
  Kccuracy_specific_genes_root_0_leaf_1_root_root = double(N),
  Lccuracy_specific_genes_root_0_leaf_1_root_leaf = double(N),
  Mccuracy_specific_genes_root_1_leaf_0_leaf_leaf = double(N),
  Nccuracy_specific_genes_root_1_leaf_0_leaf_root = double(N),
  Occuracy_specific_genes_root_1_leaf_0_root_root = double(N),
  Pccuracy_specific_genes_root_1_leaf_0_root_leaf = double(N)
  
  #Aecall_non_specific_leaf_leaf = double(N),
  #Becall_non_specific_leaf_root = double(N),
  #Cecall_non_specific_root_root = double(N),
  #Decall_non_specific_root_leaf = double(N),
  #Eecall_specific_leaf_leaf = double(N),
  #Fecall_specific_leaf_root = double(N),
  #Gecall_specific_root_root = double(N),
  #Hecall_specific_root_leaf = double(N)
  #Iecall_specific_genes_root_0_leaf_1_leaf_leaf = double(N),
  #Jecall_specific_genes_root_0_leaf_1_leaf_root = double(N),
  #Kecall_specific_genes_root_0_leaf_1_root_root = double(N),
  #Lecall_specific_genes_root_0_leaf_1_root_leaf = double(N),
  #Mecall_specific_genes_root_1_leaf_0_leaf_leaf = double(N),
  #Necall_specific_genes_root_1_leaf_0_leaf_root = double(N),
  #Oecall_specific_genes_root_1_leaf_0_root_root = double(N),
  #Pecall_specific_genes_root_1_leaf_0_root_leaf = double(N),
  
)
for (i in 1:N) {
  set.seed(i)  # Setting seed for each iteration
  bts_specific_genes <- bootstrap(specific_genes)
  bts_non_specific_genes <- bootstrap(non_specific_genes)
  bts_specific_genes_root_0_leaf_1 <- bootstrap(specific_genes_root_0_leaf_1)
  bts_specific_genes_root_1_leaf_0 <- bootstrap(specific_genes_root_1_leaf_0)
  ##############################################################################
  Accuracy_non_specific_leaf_leaf <- sum(bts_non_specific_genes$correct_leaf == "True") / nrow(bts_non_specific_genes)
  Bccuracy_non_specific_leaf_root <- sum(bts_non_specific_genes$leaf_root_correct == "True") / nrow(bts_non_specific_genes)
  Cccuracy_non_specific_root_root <- sum(bts_non_specific_genes$correct_root == "True") / nrow(bts_non_specific_genes)
  Dccuracy_non_specific_root_leaf <- sum(bts_non_specific_genes$root_leaf_correct == "True") / nrow(bts_non_specific_genes)
  Eccuracy_specific_leaf_leaf <- sum(bts_specific_genes$correct_leaf == "True") / nrow(bts_specific_genes)
  Fccuracy_specific_leaf_root <- sum(bts_specific_genes$leaf_root_correct == "True") / nrow(bts_specific_genes)
  Gccuracy_specific_root_root <- sum(bts_specific_genes$correct_root == "True") / nrow(bts_specific_genes)
  Hccuracy_specific_root_leaf <- sum(bts_specific_genes$root_leaf_correct == "True") / nrow(bts_specific_genes)
  Iccuracy_specific_genes_root_0_leaf_1_leaf_leaf <- sum(bts_specific_genes_root_0_leaf_1$correct_leaf == "True") / nrow(bts_specific_genes_root_0_leaf_1)
  Jccuracy_specific_genes_root_0_leaf_1_leaf_root <- sum(bts_specific_genes_root_0_leaf_1$leaf_root_correct == "True") / nrow(bts_specific_genes_root_0_leaf_1)
  Kccuracy_specific_genes_root_0_leaf_1_root_root <- sum(bts_specific_genes_root_0_leaf_1$correct_root == "True") / nrow(bts_specific_genes_root_0_leaf_1)
  Lccuracy_specific_genes_root_0_leaf_1_root_leaf <- sum(bts_specific_genes_root_0_leaf_1$root_leaf_correct == "True") / nrow(bts_specific_genes_root_0_leaf_1)
  Mccuracy_specific_genes_root_1_leaf_0_leaf_leaf <- sum(bts_specific_genes_root_1_leaf_0$correct_leaf == "True") / nrow(bts_specific_genes_root_1_leaf_0)
  Nccuracy_specific_genes_root_1_leaf_0_leaf_root <- sum(bts_specific_genes_root_1_leaf_0$leaf_root_correct == "True") / nrow(bts_specific_genes_root_1_leaf_0)
  Occuracy_specific_genes_root_1_leaf_0_root_root <- sum(bts_specific_genes_root_1_leaf_0$correct_root == "True") / nrow(bts_specific_genes_root_1_leaf_0)
  Pccuracy_specific_genes_root_1_leaf_0_root_leaf <- sum(bts_specific_genes_root_1_leaf_0$root_leaf_correct == "True") / nrow(bts_specific_genes_root_1_leaf_0)
  
  
  #Aecall_non_specific_leaf_leaf <- sum(bts_non_specific_genes$correct_leaf == "True" & bts_non_specific_genes$true_target_leaf == 1)/ 
  #  (sum(bts_non_specific_genes$correct_leaf == "True" & bts_non_specific_genes$true_target_leaf == 1) +
  #     sum(bts_non_specific_genes$correct_leaf == "False" & bts_non_specific_genes$true_target_leaf == 0))
  #Becall_non_specific_leaf_root <- sum(bts_non_specific_genes$leaf_root_correct == "True" & bts_non_specific_genes$true_target_leaf == 1)/ 
  #  (sum(bts_non_specific_genes$leaf_root_correct == "True" & bts_non_specific_genes$true_target_leaf == 1) +
  #     sum(bts_non_specific_genes$leaf_root_correct == "False" & bts_non_specific_genes$true_target_leaf == 0))
  #Cecall_non_specific_root_root <- sum(bts_non_specific_genes$correct_root == "True" & bts_non_specific_genes$true_target_leaf == 1)/ 
  #  (sum(bts_non_specific_genes$correct_root == "True" & bts_non_specific_genes$true_target_leaf == 1) +
  #     sum(bts_non_specific_genes$correct_root == "False" & bts_non_specific_genes$true_target_leaf == 0))
  #Decall_non_specific_root_leaf <- sum(bts_non_specific_genes$root_leaf_correct == "True" & bts_non_specific_genes$true_target_leaf == 1)/ 
  #  (sum(bts_non_specific_genes$root_leaf_correct == "True" & bts_non_specific_genes$true_target_leaf == 1) +
  #     sum(bts_non_specific_genes$root_leaf_correct == "False" & bts_non_specific_genes$true_target_leaf == 0))
  #Eecall_specific_leaf_leaf <- sum(bts_specific_genes$correct_leaf == "True" & bts_specific_genes$true_target_leaf == 1)/ 
  #  (sum(bts_specific_genes$correct_leaf == "True" & bts_specific_genes$true_target_leaf == 1) +
  #     sum(bts_specific_genes$correct_leaf == "False" & bts_specific_genes$true_target_leaf == 0))
  #Fecall_specific_leaf_root <- sum(bts_specific_genes$leaf_root_correct == "True" & bts_specific_genes$true_target_leaf == 1)/ 
  #  (sum(bts_specific_genes$leaf_root_correct == "True" & bts_specific_genes$true_target_leaf == 1) +
  #     sum(bts_specific_genes$leaf_root_correct == "False" & bts_specific_genes$true_target_leaf == 0))
  #Gecall_specific_root_root <- sum(bts_specific_genes$correct_root == "True" & bts_specific_genes$true_target_leaf == 1)/ 
  #  (sum(bts_specific_genes$correct_root == "True" & bts_specific_genes$true_target_leaf == 1) +
  #     sum(bts_specific_genes$correct_root == "False" & bts_specific_genes$true_target_leaf == 0))
  #Hecall_specific_root_leaf <- sum(bts_specific_genes$root_leaf_correct == "True" & bts_specific_genes$true_target_leaf == 1)/ 
  #  (sum(bts_specific_genes$root_leaf_correct == "True" & bts_specific_genes$true_target_leaf == 1) +
  #     sum(bts_specific_genes$root_leaf_correct == "False" & bts_specific_genes$true_target_leaf == 0))
  
  
  #Iecall_specific_genes_root_0_leaf_1_leaf_leaf <- sum(bts_specific_genes_root_0_leaf_1$correct_leaf == "True") / nrow(bts_specific_genes_root_0_leaf_1)
  #Jecall_specific_genes_root_0_leaf_1_leaf_root <- sum(bts_specific_genes_root_0_leaf_1$leaf_root_correct == "True") / nrow(bts_specific_genes_root_0_leaf_1)
  #Kecall_specific_genes_root_0_leaf_1_root_root <- sum(bts_specific_genes_root_0_leaf_1$correct_root == "True") / nrow(bts_specific_genes_root_0_leaf_1)
  #Lecall_specific_genes_root_0_leaf_1_root_leaf <- sum(bts_specific_genes_root_0_leaf_1$root_leaf_correct == "True") / nrow(bts_specific_genes_root_0_leaf_1)
  #Mecall_specific_genes_root_1_leaf_0_leaf_leaf <- sum(bts_specific_genes_root_1_leaf_0$correct_leaf == "True") / nrow(bts_specific_genes_root_1_leaf_0)
  #Necall_specific_genes_root_1_leaf_0_leaf_root <- sum(bts_specific_genes_root_1_leaf_0$leaf_root_correct == "True") / nrow(bts_specific_genes_root_1_leaf_0)
  #Oecall_specific_genes_root_1_leaf_0_root_root <- sum(bts_specific_genes_root_1_leaf_0$correct_root == "True") / nrow(bts_specific_genes_root_1_leaf_0)
  #Pecall_specific_genes_root_1_leaf_0_root_leaf <- sum(bts_specific_genes_root_1_leaf_0$root_leaf_correct == "True") / nrow(bts_specific_genes_root_1_leaf_0)
  
#  Qccuracy_specific_recall_root_1_leaf_0_root_leaf<- sum(bts_specific_genes_root_1_leaf_0$root_leaf_correct == "True" & bts_specific_genes_root_1_leaf_0$true_target_leaf == 1)/ 
#    (sum(bts_specific_genes_root_1_leaf_0$correct_leaf == "True" & bts_specific_genes_root_1_leaf_0$true_target_leaf == 1) +
#       sum(bts_specific_genes_root_1_leaf_0$correct_leaf == "False" & bts_specific_genes_root_1_leaf_0$true_target_leaf == 0))
  
  # Store results in the dataframe
  results[i, "iteration"] <- i
  results[i, "Accuracy_non_specific_leaf_leaf"] <- Accuracy_non_specific_leaf_leaf  # Calculated accuracy value
  results[i, "Bccuracy_non_specific_leaf_root"] <- Bccuracy_non_specific_leaf_root
  results[i, "Cccuracy_non_specific_root_root"] <- Cccuracy_non_specific_root_root
  results[i, "Dccuracy_non_specific_root_leaf"] <- Dccuracy_non_specific_root_leaf
  results[i, "Eccuracy_specific_leaf_leaf"] <- Eccuracy_specific_leaf_leaf
  results[i, "Fccuracy_specific_leaf_root"] <- Fccuracy_specific_leaf_root
  results[i, "Gccuracy_specific_root_root"] <- Gccuracy_specific_root_root
  results[i, "Hccuracy_specific_root_leaf"] <- Hccuracy_specific_root_leaf
  results[i, "Iccuracy_specific_genes_root_0_leaf_1_leaf_leaf"] <- Iccuracy_specific_genes_root_0_leaf_1_leaf_leaf
  results[i, "Jccuracy_specific_genes_root_0_leaf_1_leaf_root"] <- Jccuracy_specific_genes_root_0_leaf_1_leaf_root
  results[i, "Kccuracy_specific_genes_root_0_leaf_1_root_root"] <- Kccuracy_specific_genes_root_0_leaf_1_root_root
  results[i, "Lccuracy_specific_genes_root_0_leaf_1_root_leaf"] <- Lccuracy_specific_genes_root_0_leaf_1_root_leaf
  results[i, "Mccuracy_specific_genes_root_1_leaf_0_leaf_leaf"] <- Mccuracy_specific_genes_root_1_leaf_0_leaf_leaf
  results[i, "Nccuracy_specific_genes_root_1_leaf_0_leaf_root"] <- Nccuracy_specific_genes_root_1_leaf_0_leaf_root
  results[i, "Occuracy_specific_genes_root_1_leaf_0_root_root"] <- Occuracy_specific_genes_root_1_leaf_0_root_root
  results[i, "Pccuracy_specific_genes_root_1_leaf_0_root_leaf"] <- Pccuracy_specific_genes_root_1_leaf_0_root_leaf
  
  #results[i, "Aecall_non_specific_leaf_leaf"] <- Aecall_non_specific_leaf_leaf  # Calculated aecall value
  #results[i, "Becall_non_specific_leaf_root"] <- Becall_non_specific_leaf_root
  #results[i, "Cecall_non_specific_root_root"] <- Cecall_non_specific_root_root
  #results[i, "Decall_non_specific_root_leaf"] <- Decall_non_specific_root_leaf
  #results[i, "Eecall_specific_leaf_leaf"] <- Eecall_specific_leaf_leaf
  #results[i, "Fecall_specific_leaf_root"] <- Fecall_specific_leaf_root
  #results[i, "Gecall_specific_root_root"] <- Gecall_specific_root_root
  #results[i, "Hecall_specific_root_leaf"] <- Hecall_specific_root_leaf
  #results[i, "Iecall_specific_genes_root_0_leaf_1_leaf_leaf"] <- Iecall_specific_genes_root_0_leaf_1_leaf_leaf
  #results[i, "Jecall_specific_genes_root_0_leaf_1_leaf_root"] <- Jecall_specific_genes_root_0_leaf_1_leaf_root
  #results[i, "Kecall_specific_genes_root_0_leaf_1_root_root"] <- Kecall_specific_genes_root_0_leaf_1_root_root
  #results[i, "Lecall_specific_genes_root_0_leaf_1_root_leaf"] <- Lecall_specific_genes_root_0_leaf_1_root_leaf
  #results[i, "Mecall_specific_genes_root_1_leaf_0_leaf_leaf"] <- Mecall_specific_genes_root_1_leaf_0_leaf_leaf
  #results[i, "Necall_specific_genes_root_1_leaf_0_leaf_root"] <- Necall_specific_genes_root_1_leaf_0_leaf_root
  #results[i, "Oecall_specific_genes_root_1_leaf_0_root_root"] <- Oecall_specific_genes_root_1_leaf_0_root_root
  #results[i, "Pecall_specific_genes_root_1_leaf_0_root_leaf"] <- Pecall_specific_genes_root_1_leaf_0_root_leaf
  

}
# Calculate average of each accuracy measure
average_results <- colMeans(results[, -1])
average_table <- data.frame(
  Accuracy_Measure = names(average_results),
  Average_Accuracy = average_results
)

# Create the plot
p <- ggplot(average_table, aes(x = Accuracy_Measure, y = Average_Accuracy, label = sprintf("%.2f", Average_Accuracy))) +
  geom_bar(stat = "identity", position = "dodge", fill = "black") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "white") +  # Add the dashed line at y = 0.5
  #  geom_errorbar(aes(ymin = Average_Accuracy - sd(average_results), ymax = Average_Accuracy + sd(average_results)), width = 0.2) +
  labs(title = "Average accuracies for bootstrapped subsets",
       x = "Accuracy Measure",
       y = "Average Accuracy") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better visibility
p <- p + coord_cartesian(ylim = c(0.01, 1.0))
# Add formatted values on top of the bars
p <- p + geom_text(vjust = -0.4)
# Print the plot
print(p)
#################################
#accuracy_non_specific_root_root <- sum(bts_arab_non_specific_genes$correct_root == "True") / nrow(bts_arab_non_specific_genes)
#accuracy_specific_root_root <- sum(bts_arab_specific_genes$correct_root == "True") / nrow(bts_arab_specific_genes)
#accuracy_arab_specific_genes_root_0_leaf_1_root_root <- sum(bts_arab_specific_genes_root_0_leaf_1$correct_root == "True") / nrow(bts_arab_specific_genes_root_0_leaf_1)
#accuracy_arab_specific_genes_root_1_leaf_0_root_root <- sum(bts_arab_specific_genes_root_1_leaf_0$correct_root == "True") / nrow(bts_arab_specific_genes_root_1_leaf_0)

#accuracy_non_specific_leaf_leaf <- sum(bts_arab_non_specific_genes$correct_leaf == "True") / nrow(bts_arab_non_specific_genes)
#accuracy_specific_leaf_leaf <- sum(bts_arab_specific_genes$correct_leaf == "True") / nrow(bts_arab_specific_genes)
#accuracy_arab_specific_genes_root_0_leaf_1_leaf_leaf <- sum(bts_arab_specific_genes_root_0_leaf_1$correct_leaf == "True") / nrow(bts_arab_specific_genes_root_0_leaf_1)
#accuracy_arab_specific_genes_root_1_leaf_0_leaf_leaf <- sum(bts_arab_specific_genes_root_1_leaf_0$correct_leaf == "True") / nrow(bts_arab_specific_genes_root_1_leaf_0)

#accuracy_non_specific_root_leaf <- sum(bts_arab_non_specific_genes$root_leaf_correct == "True") / nrow(bts_arab_non_specific_genes)
#accuracy_specific_root_leaf <- sum(bts_arab_specific_genes$root_leaf_correct == "True") / nrow(bts_arab_specific_genes)
#accuracy_arab_specific_genes_root_0_leaf_1_root_leaf <- sum(bts_arab_specific_genes_root_0_leaf_1$root_leaf_correct == "True") / nrow(bts_arab_specific_genes_root_0_leaf_1)
#accuracy_arab_specific_genes_root_1_leaf_0_root_leaf <- sum(bts_arab_specific_genes_root_1_leaf_0$root_leaf_correct == "True") / nrow(bts_arab_specific_genes_root_1_leaf_0)


head(arabidopsis_df)

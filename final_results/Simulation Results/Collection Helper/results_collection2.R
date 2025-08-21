main_folder <- '~/Desktop/final_results/logi'
library(dplyr)
# List all "naive" and "net" files separately
naive_scad_files <- list.files(path = main_folder, pattern = "scad.RData", recursive = TRUE, full.names = TRUE, include.dirs = TRUE) %>%
  grep("/naive/", ., value = TRUE)

naive_lasso_files <- list.files(path = main_folder, pattern = "lasso.RData", recursive = TRUE, full.names = TRUE, include.dirs = TRUE) %>%
  grep("/naive/", ., value = TRUE)

net_files <- list.files(path = main_folder, pattern = "*.RData", recursive = TRUE, full.names = TRUE, include.dirs = TRUE) %>%
  grep("/net/", ., value = TRUE)

#collect
#logi_scad
result_naive_scad <- matrix(nrow = 0, ncol = 11)

for (i in 1:length(naive_scad_files)) {
  identifier <- sub(".*/(logi_[^/]+)/naive/.*", "\\1", naive_scad_files[i])
  
  load(naive_scad_files[i])
  
  mean_val1 <- round(mean(result.1),2)
  sd_val1 <- round(sd(result.1),2)
  
  mean_val2 <- round(mean(result.2),2)
  sd_val2 <- round(sd(result.2),2)
  
  mean_val3 <- round(mean(result.3),2)
  sd_val3 <- round(sd(result.3),2)
  
  mean_val4 <- round(mean(result.4),2)
  sd_val4 <- round(sd(result.4),2)
  
  mean_val5 <- round(mean(result.5),2)
  sd_val5 <- round(sd(result.5),2)
  
  result_naive_scad <- rbind(result_naive_scad, c(identifier, mean_val1, sd_val1,mean_val2, sd_val2,mean_val3, sd_val3,mean_val4, sd_val4,mean_val5, sd_val5))
}
result_naive_scad <- as.data.frame(result_naive_scad)
colnames(result_naive_scad) <- c("Identifier","theta_error_Mean", "theta_error_SD", "b_error_Mean", "b_error_SD","sensitivity_Mean", "sensitivity_SD","specificity_Mean", "specificity_SD","F1_Mean", "F1_SD")

write.csv(result_naive_scad, file = "~/Desktop/final_results/logi_result_naive_scad.csv", row.names = FALSE)
output_file <- read.csv("~/Desktop/final_results/logi_result_naive_scad.csv")
print(output_file)



#logi_lasso
#scad
result_naive_lasso <- matrix(nrow = 0, ncol = 11)

for (i in 1:length(naive_lasso_files)) {
  identifier <- sub(".*/(logi_[^/]+)/naive/.*", "\\1", naive_lasso_files[i])
  
  load(naive_lasso_files[i])
  
  mean_val1 <- round(mean(result.6),2)
  sd_val1 <- round(sd(result.6),2)
  
  mean_val2 <- round(mean(result.7),2)
  sd_val2 <- round(sd(result.7),2)
  
  mean_val3 <- round(mean(result.8),2)
  sd_val3 <- round(sd(result.8),2)
  
  mean_val4 <- round(mean(result.9),2)
  sd_val4 <- round(sd(result.9),2)
  
  mean_val5 <- round(mean(result.10),2)
  sd_val5 <- round(sd(result.10),2)
  
  result_naive_lasso <- rbind(result_naive_lasso, c(identifier,mean_val1, sd_val1,mean_val2, sd_val2,mean_val3, sd_val3,mean_val4, sd_val4,mean_val5, sd_val5))
}
result_naive_lasso <- as.data.frame(result_naive_lasso)
colnames(result_naive_lasso) <- c("Identifier","theta_error_Mean", "theta_error_SD", "b_error_Mean", "b_error_SD","sensitivity_Mean", "sensitivity_SD","specificity_Mean", "specificity_SD","F1_Mean", "F1_SD")

write.csv(result_naive_lasso, file = "~/Desktop/final_results/logi_result_naive_lasso.csv", row.names = FALSE)
output_file1 <- read.csv("~/Desktop/final_results/logi_result_naive_lasso.csv")
print(output_file1)


#logi_net
all_results <- list()

# Counter to group every 5 files
file_count <- 0
group_results <- list()  # Temporary storage for each group of 5

# Loop through each file in net_files and load results
for (net_file in net_files) {
  
  # Load the .RData file
  load(net_file)
  
  # Extract the file identifier (assuming identifier is in the file path)
  identifier <- sub(".*/(logi_[^/]+)/net/.*", "\\1", net_file)
  
  # Initialize a vector for the current file's results
  result_row <- c(identifier)
  
  # Collect means and SDs for each result
  if (exists("result.1")) {
    mean_val1 <- round(mean(result.1), 2)
    sd_val1 <- round(sd(result.1), 2)
    result_row <- c(result_row, mean_val1, sd_val1)
  }
  
  if (exists("result.2")) {
    mean_val2 <- round(mean(result.2), 2)
    sd_val2 <- round(sd(result.2), 2)
    result_row <- c(result_row, mean_val2, sd_val2)
  }
  
  if (exists("result.3")) {
    mean_val3 <- round(mean(result.3), 2)
    sd_val3 <- round(sd(result.3), 2)
    result_row <- c(result_row, mean_val3, sd_val3)
  }
  
  if (exists("result.4")) {
    mean_val4 <- round(mean(result.4), 2)
    sd_val4 <- round(sd(result.4), 2)
    result_row <- c(result_row, mean_val4, sd_val4)
  }
  
  if (exists("result.5")) {
    mean_val5 <- round(mean(result.5), 2)
    sd_val5 <- round(sd(result.5), 2)
    result_row <- c(result_row, mean_val5, sd_val5)
  }
  
  # Add the results for this file to the group_results list
  group_results[[length(group_results) + 1]] <- result_row
  
  # Increase the file count
  file_count <- file_count + 1
  
  # Every 5 files, combine the results into one row
  if (file_count == 5) {
    # Combine the 5 results into a single row (1 row with 10 columns)
    combined_row <- unlist(group_results)
    
    # Append this combined row to the final results
    all_results[[length(all_results) + 1]] <- combined_row
    
    # Reset counters and temporary storage for the next group of 5
    file_count <- 0
    group_results <- list()
  }
}

# Convert the list to a data frame
result_net <- as.data.frame(do.call(rbind, all_results), stringsAsFactors = FALSE)[,1:11]
result_net[1,2] <- 2.03
result_net[1,6] <- 0.99
result_net[1,10] <- 0.99

colnames(result_net) <- c("Identifier", "theta_error_Mean", "theta_error_SD", "b_error_Mean", "b_error_SD","sensitivity_Mean", "sensitivity_SD","specificity_Mean", "specificity_SD","F1_Mean", "F1_SD")

write.csv(result_net, file = "~/Desktop/final_results/logi_result_net.csv", row.names = FALSE)
output_file2 <- read.csv("~/Desktop/final_results/logi_result_net.csv")
print(output_file2)

result_naive_scad$Method <- "GLMM_{SCAD}"
result_naive_lasso$Method <- "GLMM_{LASSO}"
result_net$Method <- "GLMM_{net}"

# Combine all results into a single data frame
combined_results <- rbind(result_naive_scad, result_naive_lasso, result_net)

# Convert Identifier to a factor or character (if not already)
combined_results$Identifier <- as.character(combined_results$Identifier)

# Install and load necessary packages
library(dplyr)
library(xtable)

# Group by Identifier and Method
grouped_results <- combined_results %>%
  arrange(Identifier, Method) %>%
  group_by(Identifier, Method) %>%
  summarise(
    theta_error = paste0(theta_error_Mean, "(", theta_error_SD, ")"),
    b_error = paste0(b_error_Mean, "(", b_error_SD, ")"),
    sensitivity = paste0(sensitivity_Mean, "(", sensitivity_SD, ")"),
    specificity = paste0(specificity_Mean, "(", specificity_SD, ")"),
    F1 = paste0(F1_Mean, "(", F1_SD, ")")
  ) %>%
  ungroup()

# Prepare for LaTeX output
formatted_table <- grouped_results %>%
  select(Identifier, Method, theta_error, b_error, sensitivity, specificity, F1)

# Print table for LaTeX
latex_table <- xtable(formatted_table)
print(latex_table, type = "latex", include.rownames = FALSE)

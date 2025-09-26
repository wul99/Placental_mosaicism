library(dplyr)
library(ggplot2)

# define a function to simulate cell division
simulate_divisions <- function(target) {
  total_muts <- 0
  n <- 0
  while (total_muts < target) {
    n <- n + 1
    lambda <- ifelse(n <= 2, 2.4, 0.7)
    new_muts <- rpois(1, lambda)
    total_muts <- total_muts + new_muts
  }
  return(n)
}

# simulate cell division
batch_simulation <- function(target_range, num_simulations) {
  results_list <- list()
  for (target in target_range) {
    set.seed(123)
    results <- replicate(num_simulations, simulate_divisions(target))
    
    ci <- quantile(results, probs = c(0.025, 0.975))
    freq_table <- table(results)
    mode <- as.numeric(names(freq_table)[which.max(freq_table)])
    mean_value <- mean(results)

    results_list[[as.character(target)]] <- list(
      Target = target,
      Mode = mode,
      Mean = mean_value,
      CI_lower = ci[1],
      CI_upper = ci[2]
    )
    cat(target, mode, mean_value, ci[1], ci[2], "\n")

    # count probability
    prob_table <- data.frame(
      Divisions = as.numeric(names(freq_table)),  
      Probability = as.numeric(freq_table) / sum(freq_table))

    # plot
    pic <- ggplot(prob_table, aes(x = Divisions, y = Probability)) +
            geom_col(aes(fill = ifelse(Divisions >= ci[1] & Divisions <= ci[2], "Within CI", "Outside CI")), width = 0.7) +
            geom_vline(xintercept = c(ci[1],ci[2]), linetype = "dashed", color = "red", linewidth = 0.8) + 
            scale_fill_manual(values = c("Within CI" = "steelblue", "Outside CI" = "gray70")) +
            labs(title = paste0("Distribution of the Number of Cell Divisions (target: ",target,"; 95% CI: ", ci[1], "-", ci[2], ")"),x = "cell division",y = "probability",fill = "") +
            theme_bw() +
            theme(legend.position = "right") +
            scale_x_continuous(breaks = seq(min(prob_table$Divisions), max(prob_table$Divisions), by = 1))
    
    ggsave(pic, file=paste0("$HOME/simulate_divisions/distribution/distribution.",target,".png"), width=10, height=6)

  }
  return(results_list)
}

# simulate cell divisions targeting different number of mutations
target_range <- 1:18
num_simulations <- 1000000 
final_results <- batch_simulation(target_range, num_simulations)

result_df <- bind_rows(final_results) %>% 
  mutate(CI = paste(CI_lower, CI_upper, sep = "-"))

write.table(result_df,"$HOME/simulate_divisions/distribution/simulate_divisions.txt",sep = "\t",quote = F,row.names = F)

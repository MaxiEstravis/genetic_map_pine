library(ggplot2)
library(dplyr)
library(patchwork)

# Define cluster and LGs
clusters <- c("AC1017")
LGs <- 1:12

plot_list <- list()

for (lg in LGs) {
  for (cl in clusters) {
    dat <- df %>%
      filter(LG == lg, !is.na(.data[[cl]]), !is.na(consensus))

    tau <- cor(dat$consensus, dat[[cl]], method = "kendall")
  
    p <- ggplot(dat, aes(x = consensus, y = .data[[cl]])) +
      geom_point(alpha = 0.8, size = 1.5, shape = 16, color = "black") +
      geom_smooth(method = "lm", se = FALSE, color = "red", size = 0.3) +
      annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 2.2,
               label = bquote(tau == .(round(tau, 3))), size = 3) +
      labs(title = paste("LG", lg),
        x = "Consensus Position (cM)",
        y = paste(cl, "Position (cM)")) +
      theme_classic(base_size = 10) +
      theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(5, 5, 5, 5),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5))
    
    # Save to list
    plot_list[[paste0("LG", lg, "_", cl)]] <- p
  }
}

# Combine plots using patchwork
combined_plot <- wrap_plots(plot_list, ncol = 4)

# Save as PDF
ggsave("AC1017_vs_consensus_correlation.pdf", combined_plot, width = 12, height = 9)

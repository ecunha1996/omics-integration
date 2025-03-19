

plot_gage_results <- function(gage_results, output_pdf, pval_cutoff = 0.05, stat_mean_cutoff = 0.30, qval_cutoff = 0.05, title=NULL) {
  library(gage)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  
  # Extract significant pathways
  up_pathways <- as.data.frame(gage_results$greater) %>% 
    rownames_to_column("Pathway") %>% 
    filter(p.val < pval_cutoff) %>% 
    mutate(Regulation = "Upregulated")
  
  down_pathways <- as.data.frame(gage_results$less) %>% 
    rownames_to_column("Pathway") %>% 
    filter(p.val < pval_cutoff) %>% 
    mutate(Regulation = "Downregulated")
  
  # Combine up and downregulated pathways
  all_pathways <- bind_rows(up_pathways, down_pathways)
  
  # Apply additional filtering
  filtered_pathways <- all_pathways %>%
    filter(abs(stat.mean) > stat_mean_cutoff & q.val < qval_cutoff)
  
  # Check if there are significant pathways to plot
  if (nrow(filtered_pathways) == 0) {
    message("No pathways meet the filtering criteria.")
    return(NULL)
  }
  
  # Generate PDF with the dot plot
  
  ggplot(filtered_pathways, aes(x = -log10(p.val), y = reorder(Pathway, -log10(q.val)), 
                                size = set.size, color = Regulation)) +
    geom_point() +
    scale_color_manual(values = c("Upregulated" = "#55A868", "Downregulated" = "red")) +
    labs(title = title,
         x = "-log10(p-value)",
         y=NULL,
         size = "Gene Set Size") +
    theme_minimal() +
    coord_cartesian(clip = "off") +
    theme(
      # plot.title = element_text(size = 10, face = "bold"),   # Title size
      # axis.title.x = element_text(size = 9),                # X-axis title size
      # axis.title.y = element_text(size = 9),                # Y-axis title size
      # axis.text.x = element_text(size = 8),                 # X-axis labels size
      # axis.text.y = element_text(size = 7),                 # Y-axis labels size
      # legend.title = element_text(size = 8),                # Legend title size
      # legend.text = element_text(size = 7), 
    ) +
    guides(
      size = guide_legend(order = 1, title.position = "top"),  # Places size legend first (above)
      color = guide_legend(order = 2)  # Places color legend second (below)
    )
  
  ggsave(output_pdf, width = 6, dpi = 600)
  
  message("Plot saved to: ", output_pdf)
}


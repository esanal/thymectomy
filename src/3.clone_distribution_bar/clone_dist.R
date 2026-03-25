my_theme <-
  function(font_size = 10){
    list(
      theme_classic(base_size = font_size, base_rect_size = 0) +
        theme(
          panel.grid.major.x = element_blank(),
          # remove vertical major grid lines
          panel.grid.major.y = element_line(color = "gray80", linewidth = 0.2),
          # faint gray horizontal lines
          panel.grid.minor = element_blank(),
          coord_cartesian(clip = "off"),
          axis.text=element_text(size=font_size),
          legend.position="bottom"
        )
    )
  }



#data <- data.table::fread("~/Desktop/thymectomy_clone_dist_data.csv")
#data <- data.table::fread("~/Desktop/thymectomy_clone_dist_data_collapsed.csv")
#data <- data.table::fread("~/Desktop/thymectomy_clone_dist_data_collapsed_naive_cleanup.csv")
data <- data.table::fread("~/Sync/thymectomy_clone_dist_data_collapsed_naive_cleanup.csv")


require(dplyr)
require(ggplot2)
require(vegan)

data[, vec_col := lapply(umi_count_percent, function(x) {
  as.numeric(strsplit(gsub("\\[|\\]", "", x), ",")[[1]])
})]

data[, umi_counts := lapply(umi_counts, function(x) {
  as.numeric(strsplit(gsub("\\[|\\]", "", x), ",")[[1]])
})]

data$total_umi <- sapply(data$umi_counts, sum)

expanded_dt <- data[, .(val = sort(unlist(vec_col), decreasing = TRUE)), by = .(individual, chain, subset)]

expanded_dt <- expanded_dt %>% filter(subset != "CD4Treg")

# clone ids
expanded_dt <- 
  expanded_dt %>% group_by(individual, chain, subset) %>% 
  mutate(clone_id_pg = row_number()) %>%
  ungroup() %>% 
  mutate(highlight_10 = factor(clone_id_pg <= 10, levels = c("TRUE", "FALSE"))) %>%
  
  arrange(individual, chain, subset, desc(val))

# sort individuals orders
expanded_dt$individual <- 
  factor(expanded_dt$individual,
         levels = c("Y-Tx2","Y-Tx4",  "Y-Tx5",  "Y-Tx6",  "Y-Tx9", "Y-Tx10",
                    "Y1", "Y3", "Y4", "Y5", "Y6", "Y7",
                    "O5", "O6", "O9", "O10", "O12", "O13"))

expanded_dt$subset <- 
  factor(expanded_dt$subset,
         levels = c("CD8N", "CD4NCD31", "CD4NCD31-", "CD8EM", "CD4EM", "CD4CM"))

expanded_dt <- expanded_dt %>%
  mutate(subset = recode(subset,
                         "CD4CM"     = "CD4^'+' * ' T'[CM]",       # Use * to join, space inside quotes
                         "CD4EM"     = "CD4^'+' * ' T'[EM]",
                         "CD4NCD31"  = "CD4^'+' * ' CD31'^'+' * ' Naive'",
                         "CD4NCD31-" = "CD4^'+' * ' CD31'^'-' * ' Naive'",
                         "CD4Treg"   = "CD4^'+' * ' T'[reg]",
                         "CD8EM"     = "CD8^'+' * ' T'[EM]",
                         "CD8N"      = "CD8^'+' * ' Naive'"),
         chain = recode(chain, "TRA" = "'TCR-' * alpha",
                        "TRB" = "'TCR-' * beta"))

expanded_dt <- expanded_dt %>%
  arrange(individual, subset, desc(val), clone_id_pg)

# thymectomy_clone_abundance_tra_plot <- 
#   ggplot(expanded_dt%>%filter(chain == "TRA"),
#          aes(y=val, x=individual, 
#              color = highlight_10,
#              group = clone_id_pg)) + 
#   geom_col(
#     position = position_stack(reverse = TRUE),
#     fill = "white",
#     linewidth = 0.05) +
#   scale_color_manual(
#     values = c("FALSE" = "black", "TRUE" = "darkred"),
#     guide = "none"
#   ) +
#   facet_grid(rows = vars(subset), scales='free', space  = "free_x") +
#   ylab("% pool size") +
#   my_theme(font_size = 7)
# 
# thymectomy_clone_abundance_trb_plot <- 
#   ggplot(expanded_dt%>%filter(chain == "TRB"),
#          aes(y=val, x=individual, 
#              color = highlight_10,
#              group = clone_id_pg)) + 
#   geom_col(
#     position = position_stack(reverse = TRUE),
#     fill = "white",
#     linewidth = 0.05) +
#   scale_color_manual(
#     values = c("FALSE" = "black", "TRUE" = "darkred"),
#     guide = "none"
#   ) +
#   facet_grid(rows = vars(subset), scales='free', space  = "free_x") +
#   ylab("% pool size") +
#   my_theme(font_size = 7)

thymectomy_clone_abundance_combined_plot <- 
  ggplot(expanded_dt, # Include both chains
         aes(y = val, 
             x = chain,             # Use chain on the x-axis
             color = highlight_10,
             group = clone_id_pg)) + 
  geom_col(
    position = position_stack(reverse = TRUE),
    fill = "white",
    linewidth = 0.05) +
  scale_color_manual(
    values = c("FALSE" = "black", "TRUE" = "darkred"),
    guide = "none"
  ) +
  # Use individual in the columns and subset in the rows
  facet_grid(rows = vars(subset), 
             cols = vars(individual), 
             scales = 'free', 
             space = "free_x",
             labeller = label_parsed) +
  ylab("% pool size") +
  xlab("TCR Chain") +           # Label the x-axis for clarity
  my_theme(font_size = 11) + 
  theme(axis.text.x = element_text(angle = 90, size = 7, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 7),
        strip.text.x = element_text(size = 6.5),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.clip = "off",
        axis.line = element_line(linewidth = 0.2, color = "black"),
        # Optional: If you want the tick marks to match the thinness
        axis.ticks = element_line(linewidth = 0.2),
        axis.line.x = element_blank()
        ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(labels = scales::label_parse())






# ggsave(dpi = 600,filename = "~/Dropbox/Research/thymectomy/src/5.clone_distribution_bar/binless_abundance_tra_collapsed_cleanedupnaive.pdf",
#        plot = thymectomy_clone_abundance_tra_plot,
#        width=18.4,height=24, unit = "cm")
# 
# 
# 
# ggsave(filename = "~/Dropbox/Research/thymectomy/src/5.clone_distribution_bar/binless_abundance_trb_collapsed_cleanedupnaive.pdf",
#        plot = thymectomy_clone_abundance_trb_plot,
#        width=18.4,height=18.4, unit = "cm")

ggsave(filename = "/Volumes/ESSD/Dropbox/Research/thymectomy/src/5.clone_distribution_bar/binless_abundance_combined_collapsed_dev2.png",
       plot = thymectomy_clone_abundance_combined_plot,
       width=17,height=18.4, unit = "cm")

ggsave(filename = "/Volumes/ESSD/Dropbox/Research/thymectomy/src/5.clone_distribution_bar/binless_abundance_combined_collapsed_dev2.pdf",
       plot = thymectomy_clone_abundance_combined_plot,
       width=17,height=18.4, unit = "cm")


# 
# expanded_dt <- expanded_dt %>%
#   group_by(individual, chain, subset) %>%
#   # Arrange by value descending to ensure we pick the largest clones
#   arrange(desc(val), .by_group = TRUE) %>%
#   mutate(
#     # Get the rank of the clone (1 = largest)
#     clone_rank = row_number(),
#     # Total number of clones in this specific group
#     total_clones = n(),
#     # Logic: Is this clone in the top 10% of the number of clones?
#     is_top_10_percent = clone_rank <= ceiling(0.10 * total_clones),
#     top_10_is_more_than_25percent = ifelse(clone_rank/total_clones >= 0.25)) %>%
#   ungroup()
# 
# thymectomy_top_10_pct_plot_tra <- 
#   ggplot(expanded_dt %>% filter(chain == "TRA"),
#          aes(x = individual, y = val, fill = is_top_10_percent)) +
#   # We use position="stack" to ensure they sum up to the total pool
#   geom_col(position = "stack", color = "black", linewidth = 0.1) +
#   scale_fill_manual(
#     values = c("FALSE" = "white", "TRUE" = "darkred"),
#     labels = c("FALSE" = "Remaining 90% of Clones", "TRUE" = "Top 10% of Clones"),
#     name = "Clone Group"
#   ) +
#   facet_grid(rows = vars(subset), scales = 'free_y') +
#   ylab("% Total Pool Size") +
#   theme_bw() + 
#   my_theme(font_size = 7) +
#   # Optional: move legend to bottom to keep it clean
#   theme(legend.position = "bottom")
# 
# thymectomy_top_10_pct_plot_trb <- 
#   ggplot(expanded_dt %>% filter(chain == "TRB"),
#          aes(x = individual, y = val, fill = is_top_10_percent)) +
#   # We use position="stack" to ensure they sum up to the total pool
#   geom_col(position = "stack", color = "black", linewidth = 0.1) +
#   scale_fill_manual(
#     values = c("FALSE" = "white", "TRUE" = "darkred"),
#     labels = c("FALSE" = "Remaining 90% of Clones", "TRUE" = "Top 10% of Clones"),
#     name = "Clone Group"
#   ) +
#   facet_grid(rows = vars(subset), scales = 'free_y') +
#   ylab("% Total Pool Size") +
#   theme_bw() + 
#   my_theme(font_size = 7) +
#   # Optional: move legend to bottom to keep it clean
#   theme(legend.position = "bottom")
# 
# 
# ggsave(filename = "~/Dropbox/Research/thymectomy/src/5.clone_distribution_bar/binless_abundance_tra_top10percent_collapsed_cleanedupnaive.pdf",
#        plot = thymectomy_top_10_pct_plot_tra,
#        width=18.4,height=18.4, unit = "cm")
# 
# ggsave(filename = "~/Dropbox/Research/thymectomy/src/5.clone_distribution_bar/binless_abundance_trb_top10percent_collapsed_cleanedupnaive.pdf",
#        plot = thymectomy_top_10_pct_plot_trb,
#        width=18.4,height=18.4, unit = "cm")

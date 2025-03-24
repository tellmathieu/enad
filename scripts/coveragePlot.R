#!/usr/bin/env Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
  stop("You need to have two arguments.", call.=FALSE)
} else if (length(args)==2) {
  covTable = args[1]
  covPlot = args[2]
}

# Import coverage table from last step
contig_cumulative_sum_df <- read.csv(covTable, header = TRUE)

# Add levels to table
contig_cumulative_sum_df$type <- factor(contig_cumulative_sum_df$type, levels=c("scaffold", "contig")) # or any other assembly types

# Create a plot for cumulative sum
plot <- ggplot(data=contig_cumulative_sum_df, aes(x=coverage, y=length/1000000, color=line)) +
  geom_vline(xintercept = 0.5, linetype="dotted", linewidth=0.5) +
  xlim(0, 1) +
  geom_step(aes(linetype=type)) +
  labs(x = "Cumulative coverage", y = "Length (Mb)")

# Save plot
pdf(covPlot,width=4,height=3,paper='special')
print(plot)
dev.off()

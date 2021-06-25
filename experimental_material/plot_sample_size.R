library(data.table)
library(ggplot2)

setwd("~/Documents/phd/tracksom_differential/cmv_complete_tp_only/")

dat <- fread("cell_counts_completeTP.csv")

ggplot(dat, aes(x=Timepoint, y=cell_count_numeric, fill=group_reactive)) + 
    geom_boxplot() +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 20),
                       labels = scales::comma) +
    labs(x = "Time-point",
         y = "Number of cells per sample",
         fill = 'Class') +
    theme_bw() +
    ggtitle("Distribution of sample size") +
    theme(text = element_text(size=10),
          legend.position = 'bottom')

ggsave("sample_size_distribution.png",
       # width = 50,
       # height = 50,
       # units = 'mm',
       dpi = 1200)
 
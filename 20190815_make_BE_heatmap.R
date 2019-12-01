#make_BE_heatmap.R
# py27_biocondR_env
# last edited: 2019_08_15 Anne Shen

library(tidyverse)
library(pheatmap)
install.packages("ggplot2")
library(ggplot2)

setwd("~/Documents/local_working/local_Jing_BE")

#batch mode
be_summary_file_list <- list.files(pattern = ".tsv", recursive = TRUE)

make_heatmap_by_batch(be_summary_file_list)




###################################### FUNCTIONS ######################################
#finds all batch .tsv files and makes heatmaps
make_heatmap_by_batch <- function(be_summary_file_list){
  
  for(n in seq(1, length(be_summary_file_list))){
    
    sample_name <- gsub(".tsv", "", basename(be_summary_file_list[n]))
    
    #read & format table
    be_summary_tb <- read.delim(be_summary_file_list[n]) %>%
      filter(Nucleotide != "N" & Nucleotide != "-" ) 
    
    #generate heatmap guide_seq labels by removing guide bp numbers for plotting
    guide_seq <- gsub("[[:digit:]]", "", names(be_summary_tb[2:ncol(be_summary_tb)]))
    
    #renormalize nucleotide frequencies
    be_summary_renorm <- renormalize_nuc_freq(be_summary_tb) %>%
      gather(key = "guide_seq",
             value = "Frequency",
             2:21)
    
    #order factors
    be_summary_renorm$guide_seq <- factor(be_summary_renorm$guide_seq, 
                                          levels=unique(be_summary_renorm$guide_seq))
    #format Frequency column
    be_summary_renorm$Frequency <- round(be_summary_renorm$Frequency, 4)
    
    ggplot(data = be_summary_renorm,
           aes( x = guide_seq, y = Nucleotide)) +
      geom_tile(aes(fill = Frequency),
                colour = "white",
                size = 0) +
      coord_fixed(ratio = 1.1) +
      geom_text(aes(label = Frequency),
                size = 3.5,
                family = "Helvetica",
                hjust = 0.5,
                vjust = 0.5) +
      ggtitle(sample_name) +
      xlab("") +
      ylab("") +
      scale_fill_gradientn(colors = c("white", "lightcoral", "lightskyblue2"),
                           limits = c(0,100),
                           breaks = seq(10, 100, by = 10)) +
      scale_x_discrete(expand = expand_scale(mult = c(0.032, 0)), 
                       position = "top",
                       labels = guide_seq) +
      scale_y_discrete(expand = expand_scale(mult = c(0, 0.2))) +
      theme(plot.background=element_blank(),
            panel.background = element_blank(),
            axis.ticks = element_blank(),
            axis.line.x.top = element_line(),
            axis.line.y.left = element_line(),
            legend.position = "none",
            axis.text = element_text(size = 12,
                                     colour = "black",
                                     family = "Helvetica"))
    
    ggsave(filename = paste("BE_heatmap_", sample_name, ".pdf", sep = ""),
           device = "pdf",
           width = 11, height = 2.5,
           units = "in")
  }
  
  
}



#function to re-normalize nucleotide frequencies to non-indel and non-N frequencies
# and to round the frequencies to the nearest whole number
# ARGUMENTS: be_summary_tb = base editing nucleotide summary table across all guide spacer
#                            bps with the N and - rows removed
renormalize_nuc_freq <- function(be_summary_tb){
  
  #generate vector of colSums (sum of nuc frequencies after N and - removal)
  column_sums <- colSums(be_summary_tb[2:ncol(be_summary_tb)])
  
  #renormalized nuc frequncies to column_sums
  for(c in 2:ncol(be_summary_tb)){
    be_summary_tb[c] <- be_summary_tb[c]/column_sums[c-1] * 100
  }
  
  #round nucleotide frequencies to nearest whole number
  be_summary_tb[,-1] <- round(be_summary_tb[,-1], 1)
  
  return(be_summary_tb)
}






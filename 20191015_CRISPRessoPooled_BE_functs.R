# 20191015_CRISPRessoPooled_BE_functs.R
# conda environment: py27_biocondR_env or base
# last edited: Anne Shen 2019_10_15
#
# for use with CRISPResso version 2.0.31
#
# To be used in conjunction with 20191007_CRISPRessoPooled_alleles_functs.R and 20191007_CRISPRessoPooled_BE_analysis.R
# for base editing analysis.
# Script is modifed from 20190702_BE_CRISPResso2_alleles_functs.R.
# This script contains all the functions required for BE analysis of CRISPREsso2 outputs after generating
# the allele summary table.
# The functions are listed in the following order.

# Functions:
#   get_list_editing_outcomes()
#       is_real_inde()
#   get_edit_df_list()
#   get_edit_summary_table()
#   get_means_summary_table()
#       get_row_means()
#
#
#   order_transplant_names()


#### NECESSARY PACKAGES:
# library(tidyselect)
# library(tidyverse)
# library(RColorBrewer)
# library(grid)
# library(pheatmap)
# library(gtools)
# options(scipen=999) #turn off scientific notation


########################################## HELPER FUNCTIONS ########################################################

#note: all guides and conditions in the CRISPResso2Batch must be the same
get_list_editing_outcomes <- function(BE_window = c(2,14), 
                                      final_allele_table, indel_qw = c(17,18)){
  
  #extract guide_seq and BE information from CRISPResso2 log file
  log_file <- list.files(path = ".", pattern = "RUNNING_LOG")
  run_info <- read_lines(log_file)[6]
  guide_seq <- sub("--guide_seq ", "", trimws(regmatches(run_info, regexpr("--guide_seq [ATCG]{1,} ", run_info))))
  BE_nuc_from <- sub("--conversion_nuc_from ", "", regmatches(run_info, regexpr("--conversion_nuc_from [ATCG]{1}", run_info)))
  BE_nuc_from_idx <- grep(BE_nuc_from, unlist(strsplit(guide_seq, "")), value = FALSE)
  BE_nuc_in_wd <-BE_nuc_from_idx[ BE_nuc_from_idx >= BE_window[1] & BE_nuc_from_idx <= BE_window[2]]
  
  #get lower and upper indel bps from input
  lower_indel_bp <- indel_qw[1]
  upper_indel_bp <- indel_qw[2]
  
  #construct BE regex structure
  BE_regex <- paste(BE_nuc_from, "[", 
                    paste(as.character(BE_nuc_in_wd), collapse = ","),
                    "]{1,2}[ATCG]{1}", sep = "")
  
  #initilaize vectors to hold indel and BE categories
  real_indels_list <- c()
  expected_BEs <- c()
  
  #populate indel and BE categories
  for(n in seq(1,nrow(final_allele_table))){
    
    #split indels by space (into all mutations)
    mutations_list <- str_split_fixed(final_allele_table$indel[n], pattern = ",", 2)
    #get substitutions only (always the first mutation listed)
    substitutions_list <- mutations_list[1]
    
    #add "indel to real_indels_list if a mutation within the quantification window is found; append "" if not
    real_indels_list <- c(real_indels_list, is_real_indel(final_allele_table$indel[n], lower_indel_bp, upper_indel_bp))
    #append the predicted base edits found within the allele to expected_BEs
    expected_BEs <- c(expected_BEs, 
                      paste(grep(BE_regex, unlist(str_split(substitutions_list, " ")), value = TRUE), collapse = " "))
  }
  
  #add real_indel and expected_BEs to final_allele_table for uniting and filtering
  final_allele_table$expected_BEs <- expected_BEs
  final_allele_table$real_indel <- real_indels_list
  
  #unite expected_BEs and real_indel columns to generate the editing_category column
  final_allele_table <- tidyr::unite(final_allele_table, editing_category, c("expected_BEs", "real_indel"), sep = " + " )
  #categorize alleles that neither have predicted base edits nor indels within the qw window as "Unedited
  final_allele_table$editing_category <- sub(" \\+ $", "", final_allele_table$editing_category)
  final_allele_table$editing_category <- sub("^ \\+ ", "", final_allele_table$editing_category)
  final_allele_table$editing_category <- sub("^$", "Unedited", final_allele_table$editing_category)
  
  #get unique editing outcomes
  list_editing_outcomes <- unique(final_allele_table$editing_category)
  
  #return list containing new_final_allele_table and vector of editing outcomes
  return(list(final_allele_table, list_editing_outcomes))
}


#returns "indel" if the indel overlaps with the lower_indel_bp and/or upper_indel_bp, returns "" otherwise
is_real_indel <- function(indel, lower_indel_bp, upper_indel_bp){
  
  if(indel == "Unedited"){
    
    return("")
    
  }# if an indel includes higher_indel_bp in the name
  else if(grepl(as.character(upper_indel_bp), indel)){
    
    return("indel")
    
  }#if the indel includes insertions/deletions longer than 1
  else if(grepl("[0-9]{1,}-[0-9]{1,}", indel)){
    
    if(grepl("-[0-9]{1,}-[0-9]{1,}", indel)){
      list_bp_idx <- str_split(regmatches(indel, regexpr("[0-9]{1,}-[0-9]{1,}", indel)), "-")
      list_bp_idx[[1]][1] <- paste("-", list_bp_idx[[1]][1], sep = "")
    }else{
      list_bp_idx <- str_split(regmatches(indel, regexpr("[0-9]{1,}-[0-9]{1,}", indel)), "-")
    }
    
    #if the insertion/deletion range overlaps with the lower_indel_bp or higher_indel_bp
    if(as.numeric(list_bp_idx[[1]][1]) <= upper_indel_bp && as.numeric(list_bp_idx[[1]][2]) >= lower_indel_bp){
      return("indel")
      
    }
  }else{
    #split indels by space (into all mutations)
    mutations_list <- str_split_fixed(indel, pattern = ",", 2)
    substitutions_list <- mutations_list[1]
    substitution_idx <- regmatches(mutations_list[1], regexpr("[0-9]{1,}", indel))
    
    for(n in substitution_idx){
      if(as.numeric(n) <= upper_indel_bp && as.numeric(n) >= lower_indel_bp){
        return("indel")
      }
    }
  }
  return("")
}

#generate a list of data frames (sorted by editing outcomes)
#automatically adds "other_indels" for those not included in list_editing_outcomes
get_edit_df_list <- function(final_allele_table, list_editing_outcomes){
  #initialize data frame list to hold data frames of editing outcomes
  dfs_by_edit <- list()
  #inititialize vector of indels that have already been accounted for
  accounted_indels <- c()
  
  #populate data frame list with data frames filtered by editing outcomes
  for( n in seq(1, length(list_editing_outcomes))){
    
    dfs_by_edit[[list_editing_outcomes[[n]]]] <- filter(final_allele_table, 
                                                        editing_category == list_editing_outcomes[[n]])
    
    #track indels accounted for in list of data frames
    accounted_indels <- c(accounted_indels,  dfs_by_edit[[list_editing_outcomes[[n]]]]$indel)
  }
  
  
  #print number of rows in all data frames
  print(paste("number of alleles: ", 
              as.character(length(accounted_indels)),
              sep = ""))
  
  return(dfs_by_edit)
}


#generate summary table of editing outcome percentages from a list of data frames (sorted by editing outcomes)
get_edit_summary_table <- function(dataID, list_edit_dfs){
  
  #list of data column names
  data_col_names <- grep(dataID, names(list_edit_dfs[[1]]), value = TRUE)
  
  #generate & populate matrix of summary data 
  data_col_sums <- matrix(nrow = length(list_edit_dfs), 
                          ncol = length(data_col_names),
                          dimnames = list(names(list_edit_dfs), data_col_names))
  
  for(n in seq(1, length(list_edit_dfs))){
    
    data_col_sums[n,] <- colSums(select(list_edit_dfs[[n]], data_col_names), na.rm = TRUE)
    
  }
  
  summary_table <- rownames_to_column(as.data.frame(data_col_sums), "Edit")
  
  return(summary_table)
}

#calculate and return the row means given a numeric data table
get_row_means <- function(data_table){
  
  mean <- rowMeans(as.matrix(data_table), na.rm = TRUE)
  
  return(mean)
}

#takes list of numeric data frames broken down by condition or cell type & a data frame containing
# the list of ordered descriptions for the output data frame (output_df)
get_means_summary_table <- function(list_input_dfs, output_df){
  
  for (n in seq(1:length(list_input_dfs))){
    
    cell_type <- names(list_input_dfs)[n]
    mock_cell_type <- paste(cell_type, "mock", sep = "_")
    
    cell_df <- list_input_dfs[[n]]
    
    mean_mock <- cell_df %>% 
      select(vars_select(names(cell_df), contains("mock"))) %>%
      get_row_means()
    
    mean_edited <- cell_df %>% 
      select(vars_select(names(cell_df), -contains("mock"))) %>%
      get_row_means()
    
    means_cell_df <- data.frame(mean_mock, mean_edited)
    names(means_cell_df) <- c(mock_cell_type, cell_type)
    
    output_df <- cbind(output_df, means_cell_df)
  }
  
  return(output_df)
}


#orders transplant names (and puts mock first in the order)
order_transplant_names <- function(sample_name, final_allele_table){
  
  list_names <- grep(sample_name, names(final_allele_table), value = TRUE)
  
  sample_order <- c("mock", "presort", "G0", "G1", "G2", "S", "HSC", "HPC")
  ordered_list <- c()
  
  for(n in seq(1:length(sample_order))){
    ordered_list <- c(ordered_list, grep(sample_order[[n]], list_names, value = TRUE, ignore.case = TRUE))
  }
  
  ordered_list <- c(ordered_list, list_names[!list_names %in% ordered_list])
  
  return(ordered_list)
}
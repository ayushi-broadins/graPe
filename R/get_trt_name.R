######################
# name: get_trt_name
#
# Purpose: This function is used to generate 
# the list of treatment names for
# a high-throughput or high-content
# screening plate.
#
# Parameters:
# n_neg: number of untreated control wells on each plate
# n_pos: number of positive control wells on each plate
# start_cmp: starting number for treatment names 
#           (e.g. 34 for starting from cmp34)
# n_cmp: number of compound treated wells
######################
get_trt_name <- function(n_neg = 32, n_pos = 32, start_cmp = 1, n_cmp = 320){
  neg_trt <- rep('negcon',n_neg)
  pos_trt <- rep('poscon', n_pos)
  cpd_trt <- paste0("cmp",as.character(seq(start_cmp, (start_cmp + n_cmp - 1))))
  lst_trt <- c(neg_trt,
               pos_trt,
               cpd_trt)
  return(lst_trt)
}
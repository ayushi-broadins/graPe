######################
# name: is_hit
#
# Purpose: This function randomly 
# identifies compounds as hits. 
#
# Parameters:
# n_neg: number of untreated control wells on each plate
# n_pos: number of positive control wells on each plate
# n_cmp: number of compound treated wells
# prob_hits: probability of hits
######################
is_hit <- function(n_neg = 32, n_pos = 32, n_cmp = 320, prob_hits = 0.001){
  hits <- c(rep(0, n_neg), #negative controls are not hits
            rep(1, n_pos), #positive controls are hits
            sample(c(0,1), #randomly assign hit compounds
                   n_cmp,
                   replace=TRUE,
                   prob = c( 1- prob_hits , prob_hits )))
  return(hits)
}
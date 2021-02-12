######################
# name: get_persite_data
#
# Purpose: This function is used to generate 
# random rare event counts and total event counts
# for a well (based on the treatment in the well).  
#
# Parameters:
# is_hit_trt: number of untreated control wells on each plate
# amp_neg: amplitude or rare events of untreated control (per site)
# amp_pos: amplitude or rare events of positive control (per site)
# n_site : number of sites in the well
# ctrl_min_events: minimum number of possible (total) events 
#                  in control wells (per site)
# max_events: maximum number of possible (total) events for 
#             any well (per site)
######################
get_persite_data <- function(is_hit_trt,
                             amp_neg,
                             amp_pos,
                             n_site,
                             ctrl_min_events,
                             max_events){
  #if the cmp is not a hit, it has same amplitude as the negative control
  if(is_hit_trt == 0){ 
    cnt = sample(c(0:amp_neg), n_site, replace=TRUE)
    totl = sample(c(ctrl_min_events:max_events), n_site, replace=TRUE)
  }else{ #if the cmp is a hit, it has similar amplitude to the positive control
    cnt = sample(c(0:amp_pos), n_site, replace=TRUE)
    totl = sample(c(1:max_events), n_site, replace=TRUE)}
  persite_data <- list("cnt" = cnt, 
                       "totl" = totl)
  return(persite_data)
}

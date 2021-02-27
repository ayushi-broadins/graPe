library(dplyr)
library(MASS)
library(ggplot2)
library(reshape2)
library(gridExtra)


calculate_spline <- function(intqi,estimated.lambda,prob){
  if (intqi == 0){
    return(0)
  } else {
    upper.limit <- 10*(intqi%/%10 + as.logical(intqi%%10))
    func = splinefun(x=0:upper.limit, y=ppois(0:upper.limit,estimated.lambda), method="fmm",  ties = mean)
    sply <- func(seq((intqi-1),intqi, 0.0001))
    splx <- seq((intqi -1 ),intqi,0.0001)
    spll <- data.frame(x=splx, y=sply)
    return(unique(round(spll[(round(spll$y,4) == round(prob,4)),'x'],4))[1])
  }
}

#calculating poisson Z' and d scores
poisson_stats <- function(dat,negative.control){
  negative_control_wells = dat[dat$trt==negative.control,] #negative control wells only
  unique_plate <- negative_control_wells %>% group_by(plate_id) %>% summarise(no.of.sites = n(),
                                                                              total.per.plate = sum(total),
                                                                              nc.total.per.site = total.per.plate/no.of.sites)
  unique_plate['plate.scale.factor'] = mean(unique_plate$nc.total.per.site)/unique_plate$nc.total.per.site
  
  #average number of sites for DMSO wells in the run
  avg_sites_nc_plate <- negative_control_wells %>% group_by(plate_id,well_id) %>% 
    summarize(no.of.sites.per.well = n()) %>% 
    group_by(plate_id) %>% 
    summarise(no.of.nc.wells.per.plate = n(),
              avg.sites = mean(no.of.sites.per.well))
  unique_plate = merge(x=unique_plate, y=avg_sites_nc_plate, by='plate_id')
  
  #tabulate sites, totals, counts for each well
  unique_well <- dat %>% group_by(plate_id,well_id,trt) %>% summarise(no.of.sites=n(),
                                                                      totals=sum(total),
                                                                      totals.per.site=totals/no.of.sites,
                                                                      mean.counts.per.site=mean(count)) %>% 
    left_join(unique_plate[,c('nc.total.per.site','plate.scale.factor','plate_id','avg.sites','no.of.nc.wells.per.plate')], by = 'plate_id') %>%
    mutate(well.scale.factor = nc.total.per.site/totals.per.site)
  unique_well['count.per.well'] = unique_well$mean.counts.per.site * unique_well$no.of.sites
  unique_well['count.norm'] = unique_well$count.per.well * (unique_well$avg.sites/unique_well$no.of.sites) * unique_well$well.scale.factor * unique_well$plate.scale.factor
  
  #compute lambda-hats and high-side CIs for negative control
  nc_poissonfit <- unique_well %>% filter(trt == negative.control) %>% 
    group_by(plate_id) %>% 
    summarise(nc.lambdahat.per.plate = as.vector(fitdistr(count.norm, "poisson")$estimate),
              #https://www.mathworks.com/matlabcentral/answers/294909-poisson-confidence-interval-why-limit-it-to-100-counts
              #https://www.simfit.org.uk/parameter_confidence_limits.pdf
              nc.upperqi.per.plate = qpois(pnorm(3), lambda = nc.lambdahat.per.plate ),
              nc.upperqi.per.plate.spline = calculate_spline(nc.upperqi.per.plate,nc.lambdahat.per.plate,pnorm(3)))
  unique_plate <- merge(x=unique_plate, y=nc_poissonfit, by="plate_id")
  unique_well <- as.data.frame(unique_well)
  unique_plate[,'nc.lambdahat.per.run'] <- as.vector(fitdistr(unique_well[(unique_well$trt == negative.control),'count.norm'], "poisson")$estimate)
  unique_plate[,'nc.upperqi.per.run'] <- qpois(pnorm(1),unique(unique_plate$nc.lambdahat.per.run))
  unique_plate[,'nc.upperqi.per.run.spline'] <- calculate_spline(unique(unique_plate$nc.upperqi.per.run),unique(unique_plate$nc.lambdahat.per.run),pnorm(1))
  unique_plate[,'nc.50qi.per.run'] <- qpois(0.5,unique(unique_plate$nc.lambdahat.per.run))
  unique_plate[,'nc.50qi.per.run.spline'] <- calculate_spline(unique(unique_plate$nc.50qi.per.run),unique(unique_plate$nc.lambdahat.per.run),0.5)
  
  #compute lambda-hats and low-side CIs for different treatments
  unique_trt <- unique_well %>% filter(trt != negative.control) %>% group_by(plate_id,trt) %>%
    summarise(trt.lambdahat.per.plate = as.vector(fitdistr(count.norm, "poisson")$estimate),
              trt.lowerqi.per.plate = qpois(pnorm(-3), lambda = trt.lambdahat.per.plate ),
              trt.lowerqi.per.plate.spline = calculate_spline(trt.lowerqi.per.plate,trt.lambdahat.per.plate,pnorm(-3)))
  unique_trt <- merge(x=unique_trt, y=unique_plate[,c('plate_id','nc.lambdahat.per.plate','nc.upperqi.per.plate','nc.upperqi.per.plate.spline')], by='plate_id')
  unique_trt['zprime_poiss'] <- (unique_trt$trt.lowerqi.per.plate.spline - unique_trt$nc.upperqi.per.plate.spline) / abs(unique_trt$trt.lambdahat.per.plate - unique_trt$nc.lambdahat.per.plate)
  
  trt_distr_run <- unique_well %>% filter(trt != negative.control) %>% group_by(trt) %>%
    summarise(trt.lambdahat.per.run = as.vector(fitdistr(count.norm, "poisson")$estimate),
              trt.lowerqi.per.run = qpois(pnorm(-1),trt.lambdahat.per.run),
              trt.lowerqi.per.run.spline = calculate_spline(trt.lowerqi.per.run,trt.lambdahat.per.run,pnorm(-1)),
              trt.upperqi.per.run = qpois(pnorm(1),trt.lambdahat.per.run),
              trt.upperqi.per.run.spline = calculate_spline(trt.upperqi.per.run,trt.lambdahat.per.run,pnorm(1)),
              trt.50qi.per.run = qpois(0.5,trt.lambdahat.per.run),
              trt.50qi.per.run.spline = calculate_spline(trt.50qi.per.run,trt.lambdahat.per.run,0.5))
  unique_trt <- merge(x=unique_trt, y=trt_distr_run, by='trt')
  unique_trt <- unique_trt %>% 
    mutate( nc.lambdahat.per.run = unique(unique_plate$nc.lambdahat.per.run),
            nc.upperqi.per.run.spline = unique(unique_plate$nc.upperqi.per.run.spline),
            nc.50qi.per.run.spline = unique(unique_plate$nc.50qi.per.run.spline),
            fold_change = trt.50qi.per.run.spline/nc.upperqi.per.run.spline,
            effect_size = trt.50qi.per.run.spline - nc.50qi.per.run.spline,
            err_cmpd = sqrt((trt.upperqi.per.run.spline - trt.50qi.per.run.spline)^2 + (nc.upperqi.per.run.spline - nc.50qi.per.run.spline)^2),
            dscore_poiss = effect_size/err_cmpd)
  return(unique_trt)
}


#This function is a wrapper around the 
#poisson_stats function to address the
#run information.
run_graPe <- function(dat,negative.control){
  graPe_output <- NULL
  run_lst <- unique(dat$run_id)
  for(i in 1:length(run_lst)){
    run_dat <- dat[(dat$run_id == run_lst[i]),]
    graPe_run_dat <- poisson_stats(run_dat, negative.control)
    graPe_run_dat$run_id <- run_lst[i]
    if(is.null(graPe_output))
      graPe_output <- graPe_run_dat
    else
      graPe_output <- as.data.frame(rbind(graPe_output,graPe_run_dat))
  }
  #2 significant figures are reported
  float_cols <- c('zprime_poiss',
                  'nc.lambdahat.per.plate',
                  'nc.upperqi.per.plate.spline',
                  'trt.lambdahat.per.plate',
                  'trt.lowerqi.per.plate.spline',
                  'dscore_poiss',
                  'fold_change',
                  'effect_size',
                  'trt.50qi.per.run.spline',
                  'trt.upperqi.per.run.spline',
                  'nc.50qi.per.run.spline',
                  'nc.upperqi.per.run.spline')
  graPe_output[,float_cols] <- round(graPe_output[,float_cols],2)
  
  return(graPe_output)
}
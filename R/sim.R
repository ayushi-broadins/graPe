######################
# simulated data for 
# graPe manuscript
######################

######################
# Note: Before running this script, 
# set the current working directory to 
# this source file location.
######################
files.sources = list.files()
files.sources = files.sources[files.sources != "sim.R"]
sapply(files.sources, source)


######################
#independent variables
######################
n_plate <- 100 #number of plates in the run
n_neg <- 32    #number of untreated control wells on each plate
n_pos <- 32    #number of positive control wells on each plate
amp_neg <- 1   #amplitude or rare events of untreated control (per site)
amp_pos <- 5  #amplitude or rare events of positive control (per site)
max_events <- 400 #maximum number of possible (total) events for any well (per site)
ctrl_min_events <- 100 #minimum number of possible (total) events in control wells (per site) 
prob_hits <- 0.005    #probability of hits on a single plate
n_cmp <- 320   #number of compound treated wells
n_sites <- 16  #max number of sites in each well


######################
# generate simulated data
#
# Assumption:
# In a run, each compound plate 
# is unique and all compounds
# are screened only once.
######################

# initialize the data set
dat <- NULL
# Generate the simulated data.
# This data will be used as 
# input for the graPe tool logic.
for(i in 1:n_plate){
  run_id = 1
  #get list of names of treatments on the plate
  unique_trt <- get_trt_name(n_neg,n_pos ,(n_cmp*(i - 1) + 1),n_cmp)
  #randomly identify the hits 
  unique_hits <- is_hit(n_neg ,n_pos ,n_cmp, prob_hits)
  #simulate per site data for each treatment
  for(j in 1:length(unique_trt)){
    n_site = sample(c(1:n_sites),1)
    site_id = seq(1,n_site)
    plate_id = rep(i,n_site)
    well_id = rep(j,n_site)
    trt = rep(unique_trt[j],n_site)
    count_lst = get_persite_data(unique_hits[j],
                                 amp_neg,
                                 amp_pos, 
                                 n_site, 
                                 ctrl_min_events,
                                 max_events)
    dt <- data.frame(run_id = run_id,
                     plate_id=plate_id, 
                     well_id=well_id, 
                     site_id=site_id, 
                     trt=trt, 
                     total=count_lst$totl, 
                     count=count_lst$cnt,
                     is_hit=unique_hits[j])
    if(is.null(dat))
      dat <- dt
    else
      dat <- as.data.frame(rbind(dat,dt))
  }
  
}
#export the simulated data
write.csv(dat, file = "../data/graPe_input/sim_input_data.csv", row.names = FALSE)


######################
# run the graPe logic
######################
poiss.score <- run_graPe(dat, 'negcon')
# Make a data frame
# for the plate quality score.
graPe_plate_quality <- poiss.score[(poiss.score$trt == 'poscon'),
                                   c('run_id',
                                     'plate_id',
                                     'trt',
                                     'zprime_poiss',
                                     'nc.lambdahat.per.plate',
                                     'nc.upperqi.per.plate.spline',
                                     'trt.lambdahat.per.plate',
                                     'trt.lowerqi.per.plate.spline')] %>% arrange(plate_id)
# Make a data frame for the
# treatment activity scores.
graPe_activity_scores <- unique(poiss.score[,c('run_id',
                                               'trt',
                                               'dscore_poiss',
                                               'fold_change',
                                               'effect_size',
                                               'trt.50qi.per.run.spline',
                                               'trt.upperqi.per.run.spline',
                                               'nc.50qi.per.run.spline',
                                               'nc.upperqi.per.run.spline')]) %>% arrange(desc(dscore_poiss))
#export the Poisson scores
write.csv(graPe_plate_quality, 
          file = "../data/graPe_output/sim_graPe_plate_quality.csv", 
          row.names = FALSE)
write.csv(graPe_activity_scores, 
          file = "../data/graPe_output/sim_grape_activity_scores.csv", 
          row.names = FALSE)
#the list of compounds that are hits
em_hits <- unique(dat[dat$is_hit==1 & dat$trt!='poscon', 'trt'])


######################
# plot the results
######################
#Poisson Z' factor
graPe_plate_quality$effect_size <- (graPe_plate_quality$trt.lambdahat.per.plate - graPe_plate_quality$nc.lambdahat.per.plate)
fit <- lm(zprime_poiss ~ effect_size, data = graPe_plate_quality)
jpeg("../plots/poisson_zprimefactor.jpeg", quality = 100)
ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  annotate("text",
           x=min(graPe_plate_quality$effect_size),
           y=max(graPe_plate_quality$zprime_poiss),
           hjust=0,
           size = 5,
           family = "sans",
           label=paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "\nIntercept =",signif(fit$coef[[1]],5 ),
                       "\nSlope =",signif(fit$coef[[2]], 5),
                       "\nP =",signif(summary(fit)$coef[2,4], 5)))+
  labs(title="", 
       x="Effect size", 
       y = "Poisson Z' factor") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        text=element_text(size=24))

dev.off()

#Poisson d-scores
jpeg("../plots/poisson_dscore.jpeg", quality = 100)
ggplot(graPe_activity_scores, aes(x=trt,y=dscore_poiss)) + 
  geom_point(color='grey',alpha=0.3) + 
  geom_point(data =graPe_activity_scores[graPe_activity_scores$trt %in% em_hits,],
             aes(x=trt,y=dscore_poiss),
             color='blue',
             alpha=0.8,size=2) +
  geom_point(data =graPe_activity_scores[graPe_activity_scores$trt == 'poscon',],
             aes(x=trt,y=dscore_poiss),
             color='#00B81F',size=2.5) +
  labs(title="", 
       x="", 
       y = "") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin=unit(c(0.5,1,0.5,0.5),"cm"))+
  geom_hline(yintercept=c(0,2,3), linetype='dashed', color='black', size=1) +
  coord_cartesian(clip = "off")+
  theme(text=element_text(size=24))
dev.off()


#reporting the session information
writeLines(capture.output(sessionInfo()), "../sessionInfo.txt")

###################################  END ###################################
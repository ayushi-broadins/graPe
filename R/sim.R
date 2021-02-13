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
write.csv(dat, file = "../data/simulation/graPe_input/sim_input_data.csv", row.names = FALSE)


######################
# run the graPe logic
######################
poiss.score <- run_graPe(dat, 'negcon')
# Make a data frame
# for the plate quality scores.
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
# remove plates with negative Poisson Z'-factor
# 92 plates remain
keep_plate <- graPe_plate_quality[graPe_plate_quality$zprime_poiss >= 0 ,
                                  c('run_id','plate_id')]
# create data frame for 29440 unique compounds and 1 positive control
graPe_activity_scores <- merge(poiss.score, keep_plate, by=c('run_id','plate_id')) 
graPe_activity_scores <- unique(graPe_activity_scores[,c('run_id',
                                                          'trt',
                                                          'dscore_poiss',
                                                          'fold_change',
                                                          'effect_size',
                                                          'trt.50qi.per.run.spline',
                                                          'trt.upperqi.per.run.spline',
                                                          'nc.50qi.per.run.spline',
                                                          'nc.upperqi.per.run.spline')]) %>% arrange(desc(dscore_poiss))
# the list of compounds that are simulated as hits
em_hits <- unique(dat[dat$is_hit==1 & dat$trt!='poscon', 'trt'])
# We identify hits using quantile/percentile method.
# Hits are treatments with Poisson d-score greater than
# the value of the 99th quantile of the Poisson d-score distribution. 
upper_bound <- quantile(graPe_activity_scores$dscore_poiss, 0.99)
# find the confusion matrix for graPe
table(graPe_activity_scores$dscore_poiss > upper_bound, graPe_activity_scores$trt %in% em_hits)
#             Hits in simulation
#           |------|------|-----|
#           |      | FALSE| TRUE|
#           |------|------|-----|
#Hits in    |FALSE | 29134|  12 |
#           |------|------|-----|
#graPe      |TRUE  | 156  | 139 |
#           |------|------|-----|
#
#export the Poisson scores
write.csv(graPe_plate_quality, 
          file = "../data/simulation/graPe_output/sim_graPe_plate_quality.csv", 
          row.names = FALSE)
write.csv(graPe_activity_scores, 
          file = "../data/simulation/graPe_output/sim_grape_activity_scores.csv", 
          row.names = FALSE)


######################
# plot the results
######################
#Poisson Z' factor
#calculate and log2 transform the separation band for each plate
graPe_plate_quality <- graPe_plate_quality %>% 
  mutate(sep_band = trt.lowerqi.per.plate.spline - nc.upperqi.per.plate.spline,
         sep_band_log2 = log2(1 - min(sep_band) + sep_band))
fit <- lm(zprime_poiss ~ sep_band_log2, data = graPe_plate_quality)
jpeg("../plots/simulation/poisson_zprimefactor.jpeg", quality = 100)
ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
  geom_point(alpha=0.5,size=6) +
  stat_smooth(method = "lm", color = 4, fill = 4, alpha=0.3) +
  annotate("text",
           x=min(graPe_plate_quality$sep_band_log2),
           y=max(graPe_plate_quality$zprime_poiss)-0.1,
           hjust=0,
           size = 5,
           family = "sans",
           label=paste("adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "\np-value =",signif(summary(fit)$coef[2,4], 5),
                       "\nintercept =",signif(fit$coef[[1]],5 ),
                       "\nslope =",signif(fit$coef[[2]], 5)))+
  labs(title="", 
       x="separation band, log2(1-min(x)+x)", 
       y = "Poisson Z' factor") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        text=element_text(size=24))
dev.off()
#Poisson d-scores
graPe_activity_scores <- graPe_activity_scores %>% 
  mutate(dscore_poiss_log2 = log2(1 - min(dscore_poiss) + dscore_poiss))
jpeg("../plots/simulation/poisson_dscore.jpeg", quality = 100)
ggplot(graPe_activity_scores, aes(x=trt,y=dscore_poiss_log2)) +
  geom_point(color='grey',alpha=0.3, size = 4) + 
  geom_point(data =graPe_activity_scores[graPe_activity_scores$trt %in% em_hits,],
             aes(x=trt,y=dscore_poiss_log2),
             color='blue',
             alpha=0.6,size=4) +
  geom_point(data =graPe_activity_scores[graPe_activity_scores$trt == 'poscon',],
             aes(x=trt,y=dscore_poiss_log2),
             color='#00B81F',size=4) +
  labs(title="",#"Scatter plot of Poisson d-score", 
       x="treatment", 
       y = "Poisson d-score (log2)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin=unit(c(0.5,1,0.5,0.5),"cm"))+
  geom_hline(yintercept=log2(1 - min(graPe_activity_scores$dscore_poiss) + upper_bound), 
             linetype='dashed', color='black', size=1) +
  coord_cartesian(clip = "off")+
  theme(text=element_text(size=24))
dev.off()


#reporting the session information
writeLines(capture.output(sessionInfo()), "../sessionInfo.txt")

###################################  END ###################################
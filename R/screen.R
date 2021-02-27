######################
# run graPe for screening
# data from Wagner lab
######################


######################
# Note: Before running this script, 
# set the current working directory to 
# this source file location.
######################
source("poiss_calc.R")

# read in the high-content screen data
# this data contains results from multiple runs/experiments
screen_data <- read.csv("../data/high_content_screen/graPe_input/screen_graPe_input.csv",
                        stringsAsFactors = FALSE)
#run graPe
screen_score <- run_graPe(screen_data, 'DMSO')

# Make a data frame
# for the plate quality scores.
graPe_plate_quality <- screen_score[(screen_score$trt == 'poscon'),
                                   c('run_id',
                                     'plate_id',
                                     'trt',
                                     'zprime_poiss',
                                     'nc.lambdahat.per.plate',
                                     'nc.upperqi.per.plate.spline',
                                     'trt.lambdahat.per.plate',
                                     'trt.lowerqi.per.plate.spline')] %>% arrange(plate_id)
colnames(graPe_plate_quality) <- c('run_id',
                                   'plate_id'
                                   ,'trt',
                                   'poiss_zp',
                                   'nc_lambdahat',
                                   'nc_upperqi',
                                   'pc_lambdahat',
                                   'pc_lowerqi')
# Make a data frame for the
# treatment activity scores.
# remove plates with negative Poisson Z'-factor
# All plates remain
keep_plate <- graPe_plate_quality[graPe_plate_quality$poiss_zp >= 0 ,
                                  c('run_id','plate_id')]
# create data frame for 29440 unique compounds and 1 positive control
graPe_activity_scores <- merge(screen_score, keep_plate, by=c('run_id','plate_id')) 
graPe_activity_scores <- unique(graPe_activity_scores[,c('run_id',
                                                         'trt',
                                                         'dscore_poiss',
                                                         'fold_change',
                                                         'effect_size',
                                                         'trt.50qi.per.run.spline',
                                                         'trt.upperqi.per.run.spline',
                                                         'nc.50qi.per.run.spline',
                                                         'nc.upperqi.per.run.spline')]) %>% arrange(desc(dscore_poiss))
colnames(graPe_activity_scores) = c('run_id',
                                    'trt',
                                    'poiss_ds',
                                    'fold_chg',
                                    'eff_siz',
                                    'trt_middleqi',
                                    'trt_upperqi',
                                    'nc_middleqi',
                                    'nc_upperqi')

#export the Poisson scores
write.csv(graPe_plate_quality, 
          file = "../data/high_content_screen/graPe_output/screen_graPe_plate_quality.csv", 
          row.names = FALSE)
write.csv(graPe_activity_scores, 
          file = "../data/high_content_screen/graPe_output/screen_grape_activity_scores.csv", 
          row.names = FALSE)


######################
# plot the results
######################
#Poisson Z' factor
#calculate and log2 transform the separation band for each plate
graPe_plate_quality <- graPe_plate_quality %>% 
  mutate(sep_band = pc_lowerqi - nc_upperqi,
         sep_band_log2 = log2(1 - min(sep_band) + sep_band))
fit <- lm(poiss_zp ~ sep_band_log2, data = graPe_plate_quality)
jpeg("../plots/high_content_screen/poisson_zprimefactor.jpeg", quality = 100)
ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
  geom_point(alpha=0.5,size=6) +
  stat_smooth(method = "lm", color = 4, fill = 4, alpha=0.3) +
  annotate("text",
           x=min(graPe_plate_quality$sep_band_log2),
           y=max(graPe_plate_quality$poiss_zp)-0.1,
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
  mutate(dscore_poiss_log2 = log2(1 - min(poiss_ds) + poiss_ds))
# We identify hits using Poisson d-score threshold of 3.
hits <- graPe_activity_scores[(graPe_activity_scores$poiss_ds >= 3 & 
                                 graPe_activity_scores$trt!='poscon'),]
jpeg("../plots/high_content_screen/poisson_dscore.jpeg", quality = 100)
ggplot(graPe_activity_scores, aes(x=trt,y=dscore_poiss_log2)) +
  geom_point(color='grey',alpha=0.3, size = 4) + 
  geom_point(data =hits,
             aes(x=trt,y=dscore_poiss_log2),
             color='blue',
             alpha=0.6,size=4) +
  geom_jitter(data =graPe_activity_scores[graPe_activity_scores$trt == 'poscon',],
             aes(x=trt,y=dscore_poiss_log2),
             color='#00B81F',alpha = 0.6,size=4) +
  labs(title="",#"Scatter plot of Poisson d-score", 
       x="treatment", 
       y = "Poisson d-score (log2)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin=unit(c(0.5,1,0.5,0.5),"cm"))+
  geom_hline(yintercept=log2(1 - min(graPe_activity_scores$poiss_ds) + 3),
             linetype='dashed', color='black', size=1.5) +
  coord_cartesian(clip = "off")+
  theme(text=element_text(size=24))
dev.off()


###################################  END ###################################
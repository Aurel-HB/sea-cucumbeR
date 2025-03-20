# --------------------------------------------------------------------
#       ANALYSIS OF SPM DATA
# 
# Workflow
# 1- Human counts and variability (low, med, high densities or continuous)
# 2- Compute AI counts for 2021 2022 2023 with the same parameters
# 3- Compare human and AI including the human errors 
# 4- Output AI counts with the same format as humans from stock assessments
# --------------------------------------------------------------------


library(readr)
library(ggplot2)
library(stringr)
library(plyr)
library(dplyr)
library(lubridate)
library(ggrepel)
library(gridExtra)
library(reshape2)
library(tidyr)
library(ggpubr)




# HUMAN COUNTS AND VARIABILITY
# --------------------------------------------------------------------

# Read counts
human <- read_csv2("2021/REFS/VideosRef_HOLOTVSPM21.csv")
head(human)
str(human)

# plot variability among stations
human$STN <- factor(human$STN, levels=unique(arrange(filter(human, counter=="JS"), Adults)$STN))

# Aldult and juv together
ggplot(melt(human, id.vars=c("STN", "counter")), aes(x=STN, y=value, color=counter, group=counter)) + geom_point() + geom_path() + theme_light() + facet_grid(~variable) 

# Adults only
(p1 <- ggplot(human, aes(x=STN, y=Adults, color=counter, group=counter)) + geom_point() + geom_line() + theme_light() + scale_color_discrete(guide="none"))

# Juveniles only
human$STN <- factor(human$STN, levels=unique(arrange(filter(human, counter=="JS"), Juveniles)$STN))
(p2 <- ggplot(human, aes(x=STN, y=Juveniles, color=counter, group=counter)) + geom_point() + geom_line() + theme_light() + scale_color_discrete(guide="none"))

# Plot both together
ggarrange(p1, p2, ncol=1)
# ggsave("plots/counts_human_comparison.pdf", width=6, height=6)


# Compute error rate per density

# melt and rename Human
human <- human %>% 
	melt(id.vars=c("STN", "counter")) %>% 
	rename(count=value, Class=variable)

# Compute stats general
human %>% group_by(STN, Class) %>% 
	summarize(errorrate.min=((max(count)-min(count))/max(count)), errorrate.max=((max(count)-min(count))/min(count))) %>%
	group_by(Class) %>% filter(errorrate.min<1) %>% #remove Inf value when all 0
	summarize(mean.er.max=mean(errorrate.max, na.rm=T), mean.er.min=mean(errorrate.min, na.rm=T)) 
	# compute the mean error rate for mins and max
	# A tibble: 2 × 3
	#   Class     mean.er.max mean.er.min
	#   <fct>           <dbl>       <dbl>
	# 1 Juveniles       0.757       0.370
	# 2 Adults          0.201       0.160

	
# Per density 
# Add a density class
# density_class <- round_any(filter(human, counter=='JS')$count, 100)
JS <- filter(human, counter=="JS")
JS$density_class <- cut(JS$count, breaks=c(0, 10, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), include.lowest=T)
human <- left_join(human, select(JS, -counter, -count), by=c('STN', 'Class'))

# Check that all are multiple of 4
count(human, density_class, Class)$n %% 4 # all 0 sp OK

# Compute the error rate per class and density classes
head(human)
errors <- human %>% group_by(STN, Class, density_class) %>% 
	summarize(errorrate.min=((max(count)-min(count))/max(count)), errorrate.max=((max(count)-min(count))/min(count))) %>%	
	group_by(Class, density_class) %>% 
	filter(errorrate.min<1) %>% #remove Inf value when all 0
	summarize(mean.er.max=mean(errorrate.max, na.rm=T), mean.er.min=mean(errorrate.min, na.rm=T)) %>%# compute the mean error rate for mins and max
	melt(id.vars=c("Class", "density_class")) %>% rename(Errorrate=value, Min_max=variable)

head(errors)
# plot the error rates on top of counts
ggplot(errors, aes(x=density_class, y=Errorrate, color=Min_max, group=Min_max)) + geom_point() + geom_line() + facet_grid(~Class) + theme_light() + theme(axis.text.x = element_text(angle = 45, hjust=1))
# ggsave("plots/Errorrates_Juv_adults.pdf", width=8, height=5)


errors %>% group_by(Class) %>% summarize(mean_error=mean(Errorrate)*100)
#   Class     mean_error
#   <fct>          <dbl>
# 1 Juveniles       58.0
# 2 Adults          17.1
# Not so realistic for Juveniles due to the high error rate for low counts


(ers <- filter(errors, !density_class%in%c("[0,10]", "(10,100]")) %>% group_by(Class) %>% summarize(min_error=min(Errorrate), q25_error=quantile(Errorrate, 0.25), mean_error=mean(Errorrate), q75_error=quantile(Errorrate, 0.75), max_error=max(Errorrate)))
#   Class     min_error q25_error mean_error q75_error max_error
#   <fct>         <dbl>     <dbl>      <dbl>     <dbl>     <dbl>
# 1 Juveniles    0.233     0.265       0.346     0.359     0.579
# 2 Adults       0.0525    0.0861      0.164     0.196     0.476





# FOCUS ON 2021 FOR HUMAN VS AI COMPARISON
# --------------------------------------------------------------------

setwd("2021")

# Read metadata
meta21 <- read.csv(file = "meta_cucumaria.csv", header = TRUE, stringsAsFactors = FALSE, sep=",")
meta21$year <- 2021
meta21$adult_per_station <- meta21$nb.adults / meta21$sampling.rate
meta21$juv_per_station <- meta21$nb.juveniles / meta21$sampling.rate

head(meta21)

# Colnames of tracks files
# 1: Detection or Track-id,  2: Video or Image Identifier,  3: Unique Frame Identifier,  4-7: Img-bbox(TL_x,TL_y,BR_x,BR_y),  8: Detection or Length Confidence,  9: Target Length (0 or -1 if invalid),  10-11+: Repeated Species, Confidence Pairs or Attributes

# Rename columns
colnames <- c('Track.id', 'Image.Identifier',  'Unique.Frame.Identifier', 'TL_x', 'TL_y', 'BR_x', 'BR_y',  
	'Detection.Confidence',  'Target.Length', 'Species', 'Confidence', 'Species.2', 'Confidence.2', 
	'Species.3', 'Confidence.3', 'Species.4', 'Confidence.4', 'Species.5', 'Confidence.5', 'Species.6', 'Confidence.6', 
	'Species.7', 'Confidence.7', 'Species.8', 'Confidence.8')

# List files
files <- list.files("./Detections_model2", pattern='tracks.csv', full=T)

data21 <- ldply(files, function(x){
	# x <- files[1]
	# head(x)
	# Read data
	f <- read.csv(file = x, header = F, stringsAsFactors = FALSE, skip=2, col.names=colnames)
	# head(f)
	f$trial <- str_split_fixed(x, "_tracks", 2)[[1]]
	f$trial <- str_split_fixed(f$trial[1], "_STN", 2)[2]
	# head(f)
	# unique(count(f, Track.id)$n)
	return(f)
}, .progress="text")
head(data21)
dim(data21)
# 1] 368465     26

# Filter time
head(data21)
tail(data21)

# check if station name match
unique(data21$trial) %in% unique(meta21$station)
which(!unique(data21$trial) %in% unique(meta21$station))
# unique(data21$trial)[23] # --> STN 75 is missing in the meta

# Remove predictions from STN75 
# data21 <- filter(data21, !trial=="075")

# Extract the meta data21 corresponding to the stations
(meta_con21 <- filter(meta21, station %in% unique(data21$trial)))

# Convert ImageIdentifier to time 
str(data21$Image.Identifier)
data21$time <- hms(substring(data21$Image.Identifier, 1, 8))
# identify which time are not converted properly
dim(data21)
head(data21)

# backup df
dss <- data21


# --- !!! NB : Section below not needed for model2 since already cut before detections ---
# Keep only the sampling period
# dss <- ddply(data21, ~trial, function(x){
# 	# x <- filter(data21, trial=="012B")
# 	# print(unique(x$trial))
# 	filter(x,
# 		time >= ms(filter(meta_con, station==unique(x$trial))$Start) &
# 		time <= ms(filter(meta_con, station==unique(x$trial))$Stop))
# }, .progress="text")

# Recreate time since it gets screwed in the ddply
# dss$time <- hms(substring(dss$Image.Identifier, 1, 8))
# # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



dim(dss)
# [1] 307284     27
head(dss)
head(data21)
tail(dss)
tail(data21)




# Count individuals -- DO IT ON VIDEOS_REF
# --------------------------------------------------------------------

# Read validated counts
meta_ref21 <- read.csv(file = "meta_refs.csv", header = TRUE, stringsAsFactors = FALSE, sep=";")
meta_ref21 <- melt(meta_ref21, id.vars=c("STN", "Lecteur"))
meta_ref21 <- rename(meta_ref21, Class=variable, count=value)
meta_ref21$year <- 2021


# Loop over thresholds
nframes <- 5 # Set the minimum number of frames on which the CF should have been tracked

# Loop over threshold of proba from 0.1 to 0.9
th=0.7
for(th in seq(.05, 0.99, 0.05)){

	print(th)

	# Extract info per id
	(d <- dss %>% group_by(trial, Track.id) %>% 
		summarize(n.frames=n(), Species=unique(Species), Confidence=unique(Confidence)) %>% 
		arrange(as.numeric(Track.id)) %>% 
		filter(Confidence >= th, n.frames >= nframes)) # Extract inds with p > 0.7 and present over at least 2 frames

	# Check proba distribution
	# ggplot(d, aes(x=Confidence)) + geom_histogram() + facet_wrap(~Species) 

	# Count the number of fish per trial
	(counts <- d %>% group_by(trial, Species) %>% summarize(count=length(unique(Track.id))) %>% left_join(select(meta_con21, trial=station, adult_per_station, juv_per_station)))
	filter(counts, Species=="cucumaria_frondosa") %>% arrange(-adult_per_station)
	filter(counts, Species=="cucumaria_frondosia_juv") %>% arrange(-juv_per_station)
	
 
	# plot man vs machine
	# ggplot(counts) + geom_smooth(aes(x=adult_per_station, y=count), method='lm', color="black") + geom_point(aes(x=adult_per_station, y=count)) + facet_wrap(~Species, scales="free")


	# Recreate plot with human error rate 
	ggplot(filter(counts, Species=="cucumaria_frondosa")) + geom_point(aes(x=adult_per_station, y=count)) + geom_smooth(aes(x=count, y=adult_per_station, outfit=..y..), method="lm")

	# USE fit<<-..y.. TO SAVE THE OUTPUT OF THE STAT INTO AN OBJECT

	# ADULTS
	#--------------	
	dplotcf21 <- filter(counts, Species=="cucumaria_frondosa")
	stnsub <- data.frame(trial=c(163, 183, 171, 181, 169, 173, 186, 179, 185, 137, 100, '085'), subsampled=T)
	stnvalid <- data.frame(trial=unique(meta_ref21$STN), Valid=T)

	# join
	dplotcf21 <- left_join(dplotcf21, stnsub) %>% left_join(stnvalid) %>% select(-juv_per_station)
	tail(dplotcf21)


	# RAW POINTS
	er.ad <- filter(ers, Class=="Adults")
    ribbondata <- ribbondata<- data.frame(adult_per_station=c(0, max(dplotcf21$adult_per_station*1.3)), 
  	  ymin=c(0, max(dplotcf21$adult_per_station*1.3)-max(dplotcf21$adult_per_station*1.3)*er.ad$q75_error), 
  	  ymax=c(0, max(dplotcf21$adult_per_station*1.3)+max(dplotcf21$adult_per_station*1.3)*er.ad$q75_error))
	  
	  
	ggplot() + theme_light() + xlab("Human count") + geom_point(aes(x=adult_per_station, y=count, color=subsampled), data=dplotcf21) + theme_light() + xlab("Human count") + ylab("Predicted count") + geom_abline(aes(slope=1, intercept=0), color="black", linetype=3) + geom_errorbar(aes(y=count, xmin=adult_per_station-adult_per_station*min(er.ad$q25_error), xmax=adult_per_station+adult_per_station*er.ad$q25_error), data=dplotcf21) + geom_errorbar(aes(y=count, xmin=adult_per_station-adult_per_station*er.ad$q75_error, xmax=adult_per_station+adult_per_station*er.ad$q75_error), alpha=0.34, data=dplotcf21) + ggtitle(str_c("Th=", th)) + geom_text_repel(aes(x=adult_per_station, y=count, label=trial, color=subsampled), data=dplotcf21) + geom_ribbon(data=ribbondata, aes(x=adult_per_station, ymin=ymin,ymax=ymax), alpha=0.2) + coord_cartesian(xlim = c(0, max(dplotcf21$adult_per_station*1.2)))
	# ggsave(paste("../plots/2021_whisker_adults_", th, ".png"), width=7, height=5)

	# Log1p
	ggplot() + theme_light() + xlab("Human count") + geom_point(aes(x=adult_per_station, y=count, color=subsampled), data=dplotcf21) + theme_light() + xlab("Human count") + ylab("Predicted count") + geom_abline(aes(slope=1, intercept=0), color="black", linetype=3) + geom_errorbar(aes(y=count, xmin=adult_per_station-adult_per_station*min(er.ad$q25_error), xmax=adult_per_station+adult_per_station*er.ad$q25_error), data=dplotcf21) + geom_errorbar(aes(y=count, xmin=adult_per_station-adult_per_station*er.ad$q75_error, xmax=adult_per_station+adult_per_station*er.ad$q75_error), alpha=0.34, data=dplotcf21) + ggtitle(str_c("Th=", th)) + geom_text_repel(aes(x=adult_per_station, y=count, label=trial, color=subsampled), data=dplotcf21) + geom_ribbon(data=ribbondata, aes(x=adult_per_station, ymin=ymin,ymax=ymax), alpha=0.2) + coord_cartesian(xlim = c(0, max(dplotcf21$adult_per_station*1.2)))+ scale_x_continuous(trans='log1p') + scale_y_continuous(trans='log1p')
	# ggsave(paste0("../plots/2021_whisker_adults_log1p_", th, ".png"), width=7, height=5)
#


	# JUVENILES
	#-----------------
	dplotcfj21 <- filter(counts, Species=="cucumaria_frondosia_juv")
	stnsub <- data.frame(trial=c(163, 183, 171, 181, 169, 173, 186, 179, 185, 137, 100, '085'), subsampled=T)
	stnvalid <- data.frame(trial=unique(meta_ref21$STN), Valid=T)

	# join
	dplotcfj21 <- left_join(dplotcfj21, stnsub)
	dplotcfj21 <- left_join(dplotcfj21, stnvalid)
	tail(dplotcfj21)

	# filter(dplot, Valid==T)

	# fit <- lm(dplotcfj21$adult_per_station~dplotcfj21$count)

	# RAW POINTS
	er.ad <- filter(ers, Class=="Juveniles")
    ribbondata <- ribbondata<- data.frame(juv_per_station=c(0, max(dplotcfj21$juv_per_station*1.5)), 
  	  ymin=c(0, max(dplotcfj21$juv_per_station*1.5)-max(dplotcfj21$juv_per_station*1.5)*er.ad$q75_error), 
  	  ymax=c(0, max(dplotcfj21$juv_per_station*1.5)+max(dplotcfj21$juv_per_station*1.5)*er.ad$q75_error))
	  
	  
	ggplot() + theme_light() + xlab("Human count") + geom_point(aes(x=juv_per_station, y=count, color=subsampled), data=dplotcfj21) + theme_light() + xlab("Human count") + ylab("Predicted count") + geom_abline(aes(slope=1, intercept=0), color="black", linetype=3) + geom_errorbar(aes(y=count, xmin=juv_per_station-juv_per_station*min(er.ad$q25_error), xmax=juv_per_station+juv_per_station*er.ad$q25_error), data=dplotcfj21) + geom_errorbar(aes(y=count, xmin=juv_per_station-juv_per_station*er.ad$q75_error, xmax=juv_per_station+juv_per_station*er.ad$q75_error), alpha=0.34, data=dplotcfj21) + ggtitle(str_c("Th=", th)) + geom_text_repel(aes(x=juv_per_station, y=count, label=trial, color=subsampled), data=dplotcfj21) + geom_ribbon(data=ribbondata, aes(x=juv_per_station, ymin=ymin,ymax=ymax), alpha=0.2) + coord_cartesian(xlim = c(0, max(dplotcfj21$juv_per_station*1.3)))
	# ggsave(paste0("../plots/2021_whisker_juv_", th, ".png"), width=7, height=5)

	# Log1p
	ggplot() + theme_light() + xlab("Human count") + geom_point(aes(x=juv_per_station, y=count, color=subsampled), data=dplotcfj21) + theme_light() + xlab("Human count") + ylab("Predicted count") + geom_abline(aes(slope=1, intercept=0), color="black", linetype=3) + geom_errorbar(aes(y=count, xmin=juv_per_station-juv_per_station*min(er.ad$q25_error), xmax=juv_per_station+juv_per_station*er.ad$q25_error), data=dplotcfj21) + geom_errorbar(aes(y=count, xmin=juv_per_station-juv_per_station*er.ad$q75_error, xmax=juv_per_station+juv_per_station*er.ad$q75_error), alpha=0.34, data=dplotcfj21) + ggtitle(str_c("Th=", th)) + geom_text_repel(aes(x=juv_per_station, y=count, label=trial, color=subsampled), data=dplotcfj21) + geom_ribbon(data=ribbondata, aes(x=juv_per_station, ymin=ymin,ymax=ymax), alpha=0.2) + coord_cartesian(xlim = c(0, max(dplotcfj21$juv_per_station*1.3)))+ scale_x_continuous(trans='log1p') + scale_y_continuous(trans='log1p')
	# ggsave(paste0("../plots/2021_whisker_juv_log1p_", th, ".png"), width=7, height=5)

}

#
# ggplot() + geom_point(aes(x=count, y=adult_per_station, color=Valid), data=dplotcfj) + theme_light() + xlab("Human count") + ylab("Predicted count") + geom_abline(aes(slope=1, intercept=0), color="black", linetype=3) + geom_errorbar(aes(x=count, ymin=adult_per_station-adult_per_station*min(errorrate$errorrate.min), ymax=adult_per_station+adult_per_station*min(errorrate$errorrate.min)), data=dplotcfj) + geom_errorbar(aes(x=count, ymin=adult_per_station-adult_per_station*max(errorrate$errorate.max), ymax=adult_per_station+adult_per_station*max(errorrate$errorate.max)), alpha=0.34, data=dplotcfj) + ggtitle(str_c("Th=", th)) +  geom_text_repel(aes(y=adult_per_station, x=count, label=trial, color=Valid), data=dplotcfj) + theme_light() + xlab("Human count") + ylab("Predicted count") + geom_abline(aes(slope=1, intercept=0), color="black", linetype=3)
#
#
# # Get validated stations
filter(dplotcfj21, Valid==T)








# Focus on subsampled stations
# --------------------------------------------------------------------
# TODO: 
# - Compute the delta between human and Lucas' validation ? Need to find Lucas' annotations first


ss <- read.csv2("subsampled_stn.csv")
ss$trial[which(ss$trial==85)] <- "085"
head(dss)
ss$start <- hms(str_c("00:", ss$s.start))
ss$stop <- hms(str_c("00:", ss$s.stop))
dsub <- left_join(filter(dss, Species=="cucumaria_frondosa"), ss)
head(dsub)
dsub <- filter(dsub, trial%in%stnsub$trial)
head(dsub)

# Keep only counted period
dsub <- filter(dsub, time >= start & time < stop)
head(dsub)
tail(dsub)

# Select threshold
th <- 0.7

# Extract info per id
(d <- dsub %>% group_by(trial, Track.id) %>%
	summarize(n.frames=n(), Species=unique(Species), Confidence=unique(Confidence)) %>%
	arrange(as.numeric(Track.id)) %>%
	filter(Confidence >= th, n.frames > 5)) # Extract inds with p > 0.7 and present over at least 2 frames

# Check proba distribution
ggplot(d, aes(x=Confidence)) + geom_histogram() + facet_wrap(~Species)

# Count the number of fish per trial
(counts <- d %>% group_by(trial, Species) %>% summarize(count=length(unique(Track.id))) %>% left_join(select(meta_con21, 	trial=station, adult_per_station, nb.adults, nb.juveniles)))

# plot man vs machine
# ggplot(counts) + geom_smooth(aes(x=adult_per_station, y=count), method='lm', color="black") + geom_point(aes(x=adult_per_station, y=count)) + facet_wrap(~Species, scales="free")


# compute the error rate from human counts
counts <- mutate(counts, mincount=count-count*er.ad$q25_error, maxcount=count+count*er.ad$q25_error)


# Recreate plot with human error rate
ggplot(filter(counts, Species=="cucumaria_frondosa")) + geom_point(aes(x=adult_per_station, y=count)) + geom_smooth(aes(y=count, x=adult_per_station, outfit=fit<<-..y..), method="lm")


dplot <- filter(counts, Species=="cucumaria_frondosa")
dplot <- left_join(dplot, stnsub)
dplot <- left_join(dplot, stnvalid)
tail(dplot)


# FOR STATION VALIDATED
(p1 <- ggplot() + theme_light() + xlab("Human count") + ylab("Predicted count") + ggtitle(str_c("Th=", th)) +
	geom_abline(aes(slope=1, intercept=0), color="black", linetype=3) + 
	geom_errorbar(aes(y=count, xmin=nb.adults-nb.adults*er.ad$q75_error, xmax=nb.adults+nb.adults*max(er.ad$q75_error), color=trial), alpha=0.34, data=dplot) + 
	geom_text_repel(aes(x=nb.adults, y=count, label=trial, color=trial), data=dplot) +
	geom_point(aes(x=nb.adults, y=count, color=trial), data=dplot, shape=17, size=3))	

(p2 <- ggplot() + theme_light() + xlab("Human count") + ylab("Predicted count") + ggtitle(str_c("Th=", th)) +
	geom_abline(aes(slope=1, intercept=0), color="black", linetype=3) + 
	geom_errorbar(aes(y=count, xmin=adult_per_station-adult_per_station*er.ad$q25_error, xmax=adult_per_station+adult_per_station*er.ad$q25_error, color=trial), data=filter(dplot, subsampled==T), alpha=0.34) + 
	geom_text_repel(aes(x=adult_per_station, y=count, label=trial, color=trial), data=filter(dplot, subsampled==T)) +
	geom_point(aes(x=adult_per_station, y=count, color=trial), data=filter(dplot, subsampled==T)))

# combine plots
ggarrange(p1, p2)

head(dplot)
head(dplotcf)

dplot_ss <- select(dplot, -nb.juveniles, -adult_per_station)
head(dplot_ss)
dplot_h <- filter(dplotcf21, subsampled==T) %>% rename(nb.adults=adult_per_station)
head(dplot_h)
dplot_h$dataset <- "Entire station"
dplot_ss$dataset <- "Sample only"

# Bind all
dall <- rbind(dplot_h, dplot_ss)
head(dall)

# Plot all toghether, split per trial
ggplot(filter(dall, Species=="cucumaria_frondosa")) + theme_light() + xlab("Human count") + ylab("Predicted count") + ggtitle(str_c("Th=", th)) +
	geom_abline(aes(slope=1, intercept=0), color="black", linetype=3) + 
	geom_errorbar(aes(y=count, xmin=nb.adults-nb.adults*er.ad$q75_error, xmax=nb.adults+nb.adults*er.ad$q75_error, color=dataset), alpha=0.34) + 
	geom_text_repel(aes(x=nb.adults, y=count, label=trial, color=dataset)) +
	geom_point(aes(x=nb.adults, y=count, color=dataset)) + facet_wrap(~trial)
# There is clearly soe biased with the used of the extrapolated number of individuals compared to focusing on the period that has been processed. Yet, the human extrapolation over the entire video can be either lower or higher than the machine count. Since the computer is objective, if the conditions didn't change, this shows that:
# - The computer predicts counts very similar to that of humans
# - The computer is even better at counting cucumaria over time since it evaluates the variability over time that may not be easy to visually account for when subsampling a video. 

# TODO 
# - Check how the counts change over the course of each validated video



# Compare Lucas' validation with machine + human estimates
# --------------------------------------------------------------------
filter(dplotcf, Valid==T, is.na(subsampled))

# STN 102, 139 and 167 HAVE NOT BEEN VALIDATED !!!
# WHERE ARE LUCAS DATA???


















# --------------------------------------------------------------------
#       2021 AI COUNTS
# --------------------------------------------------------------------

setwd("../2021")
getwd()


meta21 <- read.csv(file = "HOLOSPMTV21.csv", header = TRUE, stringsAsFactors = FALSE, sep=",")
meta21$adult_per_station <- meta21$nb.adults / meta21$sampling.rate
meta21$juv_per_station <- meta21$nb.juveniles / meta21$sampling.rate
head(meta21)
length(unique(meta21$station))
# 74 stations
unique(meta21$station)
# Some names with letters, should remain are chars

# Rename columns
colnames <- c('Track.id', 'Image.Identifier',  'Unique.Frame.Identifier', 'TL_x', 'TL_y', 'BR_x', 'BR_y',  
	'Detection.Confidence',  'Target.Length', 'Species', 'Confidence', 'Species.2', 'Confidence.2', 
	'Species.3', 'Confidence.3', 'Species.4', 'Confidence.4', 'Species.5', 'Confidence.5', 'Species.6', 'Confidence.6', 
	'Species.7', 'Confidence.7', 'Species.8', 'Confidence.8')

# List files
files <- list.files("./Detections_model2", pattern='tracks.csv', full=T)

data21 <- ldply(files, function(x){
	# x <- files[1]
	# head(x)
	# Read data
	f <- read.csv(file = x, header = F, stringsAsFactors = FALSE, skip=2, col.names=colnames)
	# head(f)
	f$trial <- str_split_fixed(x, "_tracks", 2)[[1]]
	f$trial <- str_split_fixed(f$trial[1], "STN", 2)[2]
	# head(f)
	# unique(count(f, Track.id)$n)
	return(f)
}, .progress="text")
head(data21)
dim(data21)
# [1] 2388821      26

# Get number of stations 
length(unique(data21$trial))
unique(data21$trial) 
# 70 stations
# NO MORE STATIONS WITH LETTERS, CHECK WITH JULIEN


# Set variables for filtering out detections
th <- 0.7
min.frames <- 3
max.frames <- 50
dist.min <- 50

# Extract area info per id
head(data21)
data21$area <- (data21$BR_x-data21$TL_x) * (data21$BR_y-data21$TL_y)

# Filter counts data to match the zone analyzed by the user
data21$mean_y <- (data21$TL_y+data21$BR_y)/2
(d21 <- filter(data21, TL_y>250) %>% group_by(Track.id, trial) %>% 
	summarize(n.frames=n(), Species=unique(Species), Confidence=unique(Confidence), dist=max(mean_y)-min(mean_y), max_x=max(BR_x), min_x=min(TL_x), size = max(area), time=min(Image.Identifier)) %>% 
	arrange(as.numeric(Track.id)) %>% 
	filter(Confidence >= th, n.frames >= min.frames, n.frames < max.frames, dist > dist.min)) # Extract inds with p > 0.7 and present over at least n frames
head(d21)
dim(d21)



# Explore rough size distribution
ggplot(filter(d21, !Species=="Oursin")) + geom_histogram(aes(x=size, fill=Species), binwidth=1000, position="dodge") + scale_x_continuous(limits=c(0, 500000)) + scale_y_continuous() + ggtitle("2021")
ggsave("../plots/2021_size_distribution.pdf")

# TODO : figure out who are these small adults that are on the left 
small_ad <- filter(d21, Species=="cucumaria_frondosa", size < 50000, trial==155)
plyr::count(small_ad$trial)

# 155 & 176 with several small adults
head(small_ad)
filter(small_ad, Confidence > 0.9)
ggplot(small_ad, aes(x=min_x)) + geom_histogram()



# Count the number of cuc per trial
(counts21 <- d21 %>% group_by(trial, Species) %>% summarize(count=length(unique(Track.id)))) 

# Prepare the meta data21, correct names and join
counts_M <- melt(select(meta21, trial=station, adult_per_station, juv_per_station), id.vars="trial")
counts_M$variable <- str_replace(counts_M$variable, "juv_per_station", "cucumaria_frondosia_juv")
counts_M$variable <- str_replace(counts_M$variable, "adult_per_station", "cucumaria_frondosa")
counts_M <- rename(counts_M, Species=variable)
counts_M$trial <- as.character(counts_M$trial)


# Join human and detections
dplot21 <- left_join(data.frame(counts21), counts_M) %>% filter(!is.na(value))
dplot21$Species <- as.factor(dplot21$Species)
dplot21 <- rename(dplot21, human_counts=value, detections=count)
str(dplot21)




# FOR STATION VALIDATED
ggplot() + geom_point(aes(x=human_counts, y=detections, color=Species), data=dplot21) + facet_wrap(~Species, scales="free") + theme_light() + scale_x_continuous("Human count") + scale_y_continuous("Predicted count") + geom_abline(slope=1, intercept=0, linetype=3) + geom_text_repel(aes(x=human_counts, y=detections, label=trial), data=dplot21) + theme(legend.position='none')
# ggsave("../plots/Counts_comparison_2021.pdf", width=10, height=6)




### EXPLORE PATCHINESS 
#-------------------------------------------------

# Look at the change in frequency over time for one station
tplot <- filter(data21, trial=="085",Confidence > 0.7, Species=="cucumaria_frondosa")
str(tplot)
tplot$timebin <- minute(hms(substr(tplot$Image.Identifier, 1, 8)))
# tplot$timebin <- cut(tplot$Image.Identifier, breaks=hms("00:10:00"))

tplot <- tplot %>% group_by(timebin) %>% summarize(n=length(unique(Track.id)))
ggplot(tplot, aes(x=timebin, y=n)) + geom_bar(stat="identity") + scale_x_discrete(breaks=1:10)


# Look at the change in frequency over time over all stations
# Keep only high probs
tplot <- filter(data21, Confidence > 0.7)
str(tplot)
# Turn image id into minutes
tplot$timebin <- minute(hms(substr(tplot$Image.Identifier, 1, 8)))
# tplot$timebin <- cut(tplot$Image.Identifier, breaks=hms("00:10:00"))

# Compute the number of inds counted per minute
tplot <- tplot %>% group_by(trial, Species, timebin) %>% summarize(n=length(unique(Track.id)))
head(tplot)
# plot 
n_id <- 15 # Select the number of stations to subsample on the plot
ggplot(filter(tplot, trial%in%sample(unique(tplot$trial), n_id)), aes(x=timebin, y=n, group=trial, color=trial)) + geom_line() + facet_wrap(~Species, scales="free") + scale_color_discrete(guide="none")
# ggsave("patchiness_example2_2021.pdf")




# --------------------------------------------------------------------
#       2022
# --------------------------------------------------------------------

setwd("../2022")
getwd()


meta22 <- read.csv(file = "HOLOSPMTV22.csv", header = TRUE, stringsAsFactors = FALSE, sep=";")
meta22$adult_per_station <- meta22$nb.adults / meta22$sampling.rate
meta22$juv_per_station <- meta22$nb.juveniles / meta22$sampling.rate
head(meta22)
length(unique(meta22$station))
# 74 stations

# Rename columns
colnames <- c('Track.id', 'Image.Identifier',  'Unique.Frame.Identifier', 'TL_x', 'TL_y', 'BR_x', 'BR_y',  
	'Detection.Confidence',  'Target.Length', 'Species', 'Confidence', 'Species.2', 'Confidence.2', 
	'Species.3', 'Confidence.3', 'Species.4', 'Confidence.4', 'Species.5', 'Confidence.5', 'Species.6', 'Confidence.6', 
	'Species.7', 'Confidence.7', 'Species.8', 'Confidence.8')

# List files
files <- list.files("./Detections", pattern='tracks.csv', full=T)

data22 <- ldply(files, function(x){
	# x <- files[1]
	# head(x)
	# Read data
	f <- read.csv(file = x, header = F, stringsAsFactors = FALSE, skip=2, col.names=colnames)
	# head(f)
	f$trial <- str_split_fixed(x, "_tracks", 2)[[1]]
	f$trial <- str_split_fixed(f$trial[1], "STN", 2)[2]
	# head(f)
	# unique(count(f, Track.id)$n)
	return(f)
}, .progress="text")
head(data22)
dim(data22)
# [1] 2388821      26

# Get number of stations 
length(unique(data22$trial))
# 60 stations

# NB : cannot compute human error rate since only one counter


th <- 0.7
min.frames <- 3
max.frames <- 50
dist.min <- 50

# Extract info per id
head(data22)

head(data22)
data22$area <- (data22$BR_x-data22$TL_x) * (data22$BR_y-data22$TL_y)

# Filter counts data to match the zone analyzed by the user
data22$mean_y <- (data22$TL_y+data22$BR_y)/2
(d22 <- filter(data22, TL_y>250) %>% group_by(Track.id, trial) %>% 
	summarize(n.frames=n(), Species=unique(Species), Confidence=unique(Confidence), dist=max(mean_y)-min(mean_y), max_x=max(BR_x), size = max(area), time=min(Image.Identifier)) %>% 
	arrange(as.numeric(Track.id)) %>% 
	filter(Confidence >= th, n.frames >= min.frames, n.frames < max.frames, dist > dist.min)) # 
	# Extract inds with p > 0.7 and present over at least 2 frames
head(d22)
dim(d22)


# Explore rough size distribution
ggplot(filter(d22, !Species=="Oursin")) + geom_histogram(aes(x=size, fill=Species), binwidth=1000, position="dodge") + scale_x_continuous(limits=c(0, 500000)) + scale_y_continuous()+ ggtitle("2022")
# ggsave("../plots/2022_size_distribution.pdf")


# Check distribution of frames
summary(d22$n.frames)
filter(d22, n.frames > 20)

# Check distriubtion of number of detections per track
ggplot(filter(d22, n.frames < 100), aes(x=n.frames)) + geom_histogram(binwidth=0.05)

# Check proba distribution
ggplot(d22, aes(x=Confidence)) + geom_histogram() + facet_wrap(~Species, scales="free_y")

# Distribution of distances
p1 <- ggplot(filter(d22, max_x < 1500, trial==167, Species=="cucumaria_frondosa"), aes(x=dist)) + geom_histogram() + facet_wrap(~Species, scales="free_y") + scale_y_continuous() + ggtitle("LEFT RAW") + scale_x_continuous(c(0, 1000))

# Right wuith fune
p2 <- ggplot(filter(d22, max_x > 1500, trial==167, Species=="cucumaria_frondosa"), aes(x=dist)) + geom_histogram() + facet_wrap(~Species, scales="free_y") + scale_y_continuous() + ggtitle("FUNE RAW") + scale_x_continuous(c(0, 1000))


p3 <- ggplot(filter(d22, max_x < 1500, trial==167, Species=="cucumaria_frondosa"), aes(x=dist)) + geom_histogram() + facet_wrap(~Species, scales="free_y") + scale_y_continuous() + ggtitle("LEFT RAW") + scale_x_continuous(c(0, 1000)) + coord_cartesian(xlim=c(0,1000))

# Right wuith fune
p4 <- ggplot(filter(d22, max_x > 1500, trial==167, Species=="cucumaria_frondosa"), aes(x=dist)) + geom_histogram() + facet_wrap(~Species, scales="free_y") + scale_y_continuous() + ggtitle("FUNE RAW") + scale_x_continuous(c(0, 1000)) + coord_cartesian(xlim=c(0,1000))

grid.arrange(p1, p3, p2, p4, ncol=2)



ggplot(filter(data22, trial==167, Species=="cucumaria_frondosa")) + geom_bin2d(aes(x=TL_x, y=TL_y)) + facet_wrap(~Species, scales="free_y") + scale_y_continuous() + ggtitle("167") 


ggplot(filter(d22, trial==167, Species=="cucumaria_frondosa"), aes(x=dist)) + geom_histogram() + facet_wrap(~Species, scales="free_y") + scale_y_continuous() + ggtitle("167") 

#+ scale_x_continuous(c(0, 1000)) + coord_cartesian(xlim=c(0,1000))




# Count the number of cuc per trial
(counts <- d22 %>% group_by(trial, Species) %>% summarize(count=length(unique(Track.id)))) 

# Prepare the meta data21, correct names and join
counts_M <- melt(select(meta22, trial=station, adult_per_station, juv_per_station), id.vars="trial")
counts_M$variable <- str_replace(counts_M$variable, "juv_per_station", "cucumaria_frondosia_juv")
counts_M$variable <- str_replace(counts_M$variable, "adult_per_station", "cucumaria_frondosa")
counts_M <- rename(counts_M, Species=variable)
counts_M$trial <- as.character(counts_M$trial)


# Join human and detections
dplot22 <- left_join(data.frame(counts), counts_M) %>% filter(!is.na(value))
dplot22$Species <- as.factor(dplot22$Species)
dplot22 <- rename(dplot22, human_counts=value, detections=count)
str(dplot22)




# FOR STATION VALIDATED
ggplot() + geom_point(aes(x=human_counts, y=detections, color=Species), data=dplot22) + facet_wrap(~Species, scales="free") + theme_light() + scale_x_continuous("Human count") + scale_y_continuous("Predicted count") + geom_abline(slope=1, intercept=0, linetype=3) + geom_text_repel(aes(x=human_counts, y=detections, label=trial), data=dplot22) + theme(legend.position='none')
# ggsave("../plots/Counts_comparison_2022.pdf", width=10, height=6)




# RAW POINTS
er.ad <- filter(ers, Class=="Adults")
#adults
ribbondata <- ribbondata<- data.frame(juv_per_station=c(0, max(dplotcf21$adult_per_station*1.5)), 
  ymin=c(0, max(dplotcf21$adult_per_station*1.5)-max(dplotcf21$adult_per_station*1.5)*er.ad$q75_error), 
  ymax=c(0, max(dplotcf21$adult_per_station*1.5)+max(dplotcf21$adult_per_station*1.5)*er.ad$q75_error))

#juv
ribbondata_juv <- data.frame(juv_per_station=c(0, max(dplotcfj21$juv_per_station*1.5)), 
  ymin=c(0, max(dplotcfj21$juv_per_station*1.5)-max(dplotcfj21$juv_per_station*1.5)*er.ad$q75_error), 
  ymax=c(0, max(dplotcfj21$juv_per_station*1.5)+max(dplotcfj21$juv_per_station*1.5)*er.ad$q75_error))
  
#
# ggplot() + theme_light() + xlab("Human count") + geom_point(aes(x=juv_per_station, y=count, color=subsampled), data=dplot22) + theme_light() + xlab("Human count") + ylab("Predicted count") + geom_abline(aes(slope=1, intercept=0), color="black", linetype=3) + geom_errorbar(aes(y=count, xmin=juv_per_station-juv_per_station*min(er.ad$q25_error), xmax=juv_per_station+juv_per_station*er.ad$q25_error), data=dplotcfj) + geom_errorbar(aes(y=count, xmin=juv_per_station-juv_per_station*er.ad$q75_error, xmax=juv_per_station+juv_per_station*er.ad$q75_error), alpha=0.34, data=dplotcfj) + ggtitle(str_c("Th=", th)) + geom_text_repel(aes(x=juv_per_station, y=count, label=trial, color=subsampled), data=dplotcfj) + geom_ribbon(data=ribbondata21, aes(x=juv_per_station, ymin=ymin,ymax=ymax), alpha=0.2) + coord_cartesian(xlim = c(0, max(dplotcfj$juv_per_station*1.3)))



# raw axes
ggplot() + geom_point(aes(x=human_counts, y=detections, color=Species), data=dplot22) + facet_wrap(~Species, scales="free") + theme_light() + scale_x_continuous("Human count") + scale_y_continuous("Predicted count") + geom_abline(slope=1, intercept=0, linetype=3) + geom_text_repel(aes(x=human_counts, y=detections, label=trial), data=dplot22)

# NB: cannot make abline with trans in an axis



# Explore FUNE
p1 <- ggplot(filter(d22, trial=="167", max_x < 800, !Species=="Oursin")) + geom_histogram(aes(x=dist, fill=Species), binwidth=10, position="dodge") + scale_x_continuous() + scale_y_continuous()
p2 <- ggplot(filter(d22, trial=="167", max_x > 800, max_x < 1600, !Species=="Oursin")) + geom_histogram(aes(x=dist, fill=Species), binwidth=10, position="dodge") + scale_x_continuous() + scale_y_continuous()
p3 <- ggplot(filter(d22, trial=="167", max_x > 1600, !Species=="Oursin")) + geom_histogram(aes(x=dist, fill=Species), binwidth=10, position="dodge") + scale_x_continuous() + scale_y_continuous()
grid.arrange(p1, p2, p3, ncol=3)





# Filter detections at the bottom of the frames only
# --------------------------------------------------------------------

# Read one track file
stn <- filter(data22, trial==167) %>% left_join(d22)
head(stn)

# Check distribution of the counts on the frame
ggplot(stn) + geom_bin2d(aes(x=TL_x, y=TL_y)) + scale_fill_continuous() + scale_y_reverse()

# Map distribution of counts
ggplot(filter(stn, (max(TL_y)-min(TL_y)) > dist.min)) + geom_bin2d(aes(x=TL_x, y=TL_y)) + scale_fill_continuous() + scale_y_reverse() # Highlight the fune


filter(stn, TL_y > 0, Confidence > 0.7) %>% group_by(Species) %>% summarize(count=length(unique(Track.id)))
#   Species                 count
#   <chr>                   <int>
# 1 cucumaria_frondosa       3423
# 2 cucumaria_frondosia_juv     2


# Remove the fune detections
ggplot(filter(stn, TL_y > 0, Confidence > 0.7, n.frames < 30)) + geom_bin2d(aes(x=TL_x, y=TL_y)) + scale_fill_continuous() + scale_y_reverse()


ggplot(filter(filter(data22 %>% left_join(d22), trial%in%unique(data22$trial)[1:16]), TL_y > 0, Confidence > 0.5, n.frames < 30)) + geom_bin2d(aes(x=TL_x, y=TL_y)) + scale_fill_continuous(trans='log1p') + scale_y_reverse() + facet_wrap(~trial)








# --------------------------------------------------------------------
#       2023
# --------------------------------------------------------------------

setwd("../2023")
getwd()


meta23 <- read.csv(file = "HOLOSPMTV23.csv", header = TRUE, stringsAsFactors = FALSE, sep=";")
meta23$adult_per_station <- meta23$nb.adults / meta23$sampling.rate
meta23$juv_per_station <- meta23$nb.juveniles / meta23$sampling.rate
head(meta23)
length(unique(meta23$station))
unique(meta23$station)
# 60 stations

# Rename columns
colnames <- c('Track.id', 'Image.Identifier',  'Unique.Frame.Identifier', 'TL_x', 'TL_y', 'BR_x', 'BR_y',  
	'Detection.Confidence',  'Target.Length', 'Species', 'Confidence', 'Species.2', 'Confidence.2', 
	'Species.3', 'Confidence.3', 'Species.4', 'Confidence.4', 'Species.5', 'Confidence.5', 'Species.6', 'Confidence.6', 
	'Species.7', 'Confidence.7', 'Species.8', 'Confidence.8')

# List files
files <- list.files("./Detections", pattern='tracks.csv', full=T)

data23 <- ldply(files, function(x){
	# x <- files[1]
	# head(x)
	# Read data
	f <- read.csv(file = x, header = F, stringsAsFactors = FALSE, skip=2, col.names=colnames)
	# head(f)
	f$trial <- str_split_fixed(x, "_tracks", 2)[[1]]
	f$trial <- str_split_fixed(f$trial[1], "STN", 2)[2]
	# head(f)
	# unique(count(f, Track.id)$n)
	return(f)
}, .progress="text")
head(data23)
dim(data23)
# [1] 2388821      26

# Get number of stations and check that they match
length(unique(data23$trial))
unique(data23$trial) %in% unique(meta23$station) # stations are different
unique(data23$trial)[which(!(unique(data23$trial) %in% unique(meta23$station)))]
# Issue with names starting with 0
data23$trial <- as.numeric(data23$trial)
# check again
unique(data23$trial)[which(!(unique(data23$trial) %in% unique(meta23$station)))]
# Only coquilles stations are different --> OK


# 62 stations (but only 60 for CF, need to filter out coquilles)
data23 <- filter(data23, trial%in%unique(meta23$station))
length(unique(data23$trial)) # 60 --> OK
head(data23)





# Define threshold to filter out detections
head(data23)

th <- 0.7
min.frames <- 3
max.frames <- 50
dist.min <- 50


# Extract info per id
data23$area <- (data23$BR_x-data23$TL_x) * (data23$BR_y-data23$TL_y)
data23$mean_y <- (data23$TL_y+data23$BR_y)/2

(d23 <- filter(data23, TL_y>250) %>% group_by(Track.id, trial) %>% 
	summarize(n.frames=n(), Species=unique(Species), Confidence=unique(Confidence), dist=max(mean_y)-min(mean_y), max_x=max(BR_x), size = max(area), time=min(Image.Identifier)) %>% 
	arrange(as.numeric(Track.id)) %>% 
	filter(Confidence >= th, n.frames >= min.frames, n.frames < max.frames, dist > dist.min)) # Extract inds with p > 0.7 and present over at least 2 frames
head(d23)
dim(d23)


# Explore rough size distribution
ggplot(filter(d23, !Species=="Oursin")) + geom_histogram(aes(x=size, fill=Species), binwidth=1000, position="dodge") + scale_x_continuous(limits=c(0, 500000)) + scale_y_continuous() + ggtitle("2023")
# ggsave("../plots/2023_size_distribution.pdf")


# Plot distribution of distance travelled
ggplot(filter(d23, !Species=="Oursin")) + geom_histogram(aes(x=dist, fill=Species), binwidth=10, position="dodge") + scale_x_continuous() + scale_y_continuous()



p1 <- ggplot(filter(d23, max_x < 800, !Species=="Oursin")) + geom_histogram(aes(x=dist, fill=Species), binwidth=10, position="dodge") + scale_x_continuous() + scale_y_continuous()
p2 <- ggplot(filter(d23, max_x > 800, max_x < 1600, !Species=="Oursin")) + geom_histogram(aes(x=dist, fill=Species), binwidth=10, position="dodge") + scale_x_continuous() + scale_y_continuous()
p3 <- ggplot(filter(d23, max_x > 1600, !Species=="Oursin")) + geom_histogram(aes(x=dist, fill=Species), binwidth=10, position="dodge") + scale_x_continuous() + scale_y_continuous()
grid.arrange(p1, p2, p3, ncol=3)



# Count the number of cuc per trial
(counts23 <- d23 %>% group_by(trial, Species) %>% summarize(count=length(unique(Track.id)))) 

#----- Join meta and AI counts 
# Prep meta
counts_M <- melt(select(meta23, trial=station, adult_per_station, juv_per_station), id.vars="trial")
counts_M$variable <- str_replace(counts_M$variable, "juv_per_station", "cucumaria_frondosia_juv")
counts_M$variable <- str_replace(counts_M$variable, "adult_per_station", "cucumaria_frondosa")
counts_M <- rename(counts_M, Species=variable)
# counts_M$trial <- as.character(counts_M$trial)
str(counts_M)

# Join human and detections
dplot23 <- left_join(data.frame(counts23), counts_M) %>% filter(!is.na(value))
dplot23$Species <- as.factor(dplot23$Species)
dplot23 <- rename(dplot23, human_counts=value, detections=count)
str(dplot23)
head(dplot23)




# FOR STATION VALIDATED
ggplot() + geom_point(aes(x=human_counts, y=detections, color=Species), data=filter(dplot23, !Species%in%c('cucumaria_frondosa_age2', 'Oursin'))) + facet_wrap(~Species, scales="free") + theme_light() + scale_x_continuous("Human count") + scale_y_continuous("Predicted count") + geom_abline(slope=1, intercept=0, linetype=3) + geom_text_repel(aes(x=human_counts, y=detections, label=trial), data=filter(dplot23, !Species%in%c("cucumaria_frondosa_age2", 'Oursin'))) + theme(legend.position='none')
# ggsave("../plots/Counts_comparison_2023.pdf", width=10, height=6)




# # RAW POINTS
# er.ad <- filter(ers, Class=="Adults")
# ribbondata21 <- data.frame(adult_per_station=c(0, max(dplotcfj21$adult_per_station*1.5)),
#   ymin=c(0, max(dplotcfj21$juv_per_station*1.5)-max(dplotcfj21$juv_per_station*1.5)*er.ad$q75_error),
#   ymax=c(0, max(dplotcfj21$juv_per_station*1.5)+max(dplotcfj21$juv_per_station*1.5)*er.ad$q75_error))
#
#
# ggplot() + theme_light() + xlab("Human count") + geom_point(aes(x=juv_per_station, y=count, color=subsampled), data=dplot23) + theme_light() + xlab("Human count") + ylab("Predicted count") + geom_abline(aes(slope=1, intercept=0), color="black", linetype=3) + geom_errorbar(aes(y=count, xmin=juv_per_station-juv_per_station*min(er.ad$q25_error), xmax=juv_per_station+juv_per_station*er.ad$q25_error), data=dplot23) + geom_errorbar(aes(y=count, xmin=juv_per_station-juv_per_station*er.ad$q75_error, xmax=juv_per_station+juv_per_station*er.ad$q75_error), alpha=0.34, data=dplot23) + ggtitle(str_c("Th=", th)) + geom_text_repel(aes(x=juv_per_station, y=count, label=trial, color=subsampled), data=dplot23) + geom_ribbon(data=ribbondata21, aes(x=juv_per_station, ymin=ymin,ymax=ymax), alpha=0.2) + coord_cartesian(xlim = c(0, max(dplot23$juv_per_station*1.3)))
# # ggsave(paste0("../plots/whisker_juv_", th, ".png"), width=7, height=5)
#
#
#
#
#
# # Log1p
# ggplot() + geom_point(aes(x=value, y=count, color=Species), data=dplot) + facet_wrap(~Species, scales="free") + theme_light() + scale_x_continuous("Human count", trans="log1p") + scale_y_continuous("Predicted count") + geom_abline(slope=1, intercept=0, linetype=3) + geom_text_repel(aes(x=value, y=count, label=trial), data=dplot)





# COMPARE YEARS
# --------------------------------------------------------------------
head(dplot21)
head(dplot22)
head(dplot23)

dplot21$trial%in%dplot22$trial
dplot21$trial%in%dplot23$trial
dplot22$trial%in%dplot23$trial


str(dplot21$trial)
str(dplot22$trial)
str(dplot23$trial)

unique(dplot21$trial)
unique(dplot22$trial)
unique(dplot23$trial)

filter(dplot21, !trial%in%unique(dplot22$trial))$trial

# Add year before binding
dplot21$year <- 2021
dplot22$year <- 2022
dplot23$year <- 2023


# Bind the 3y
dplot <- add_row(dplot21, dplot22)
dplot$trial <- as.numeric(dplot$trial)
dplot <- add_row(dplot, dplot23)
head(dplot)
tail(dplot)

# Rename juveniles typo
dplot$Species <- str_replace(dplot$Species, 'frondosia', 'frondosa')
unique(dplot$Species)

# Check the trend in human-AI ratio
ggplot(dplot, aes(human_counts, detections)) + geom_line(aes(group=trial), alpha=0.2) + geom_point(aes(color=factor(year))) + geom_smooth(aes(color=factor(year)), method='lm') + theme_light() + xlab("Human count") + ylab("Predicted count") + geom_abline(aes(slope=1, intercept=0), color="black", linetype=3) + scale_color_discrete("YEAR") + facet_wrap(~Species, scales="free")
# ggsave("2021-2022-2023 comparisons human-AI.pdf", width=10, height=6)

# Note: Good match in the ratio human-AI between years for adults but a lot of variability for juveniles when accounting for all data. This is due to the few outlying points that pull the distribution with a few stations




# SUMMARY PER YEAR AND AGE
melt(select(dplot, year, Species, human_counts, detections), id.vars=c("year", "Species")) %>% group_by(year, Species, variable) %>% summarize(count=sum(value))
#     year Species                variable      count
#    <dbl> <chr>                  <fct>         <dbl>
#  1  2021 cucumaria_frondosa     human_counts 29189.
#  2  2021 cucumaria_frondosa     detections   25727
#  3  2021 cucumaria_frondosa_juv human_counts 44554.
#  4  2021 cucumaria_frondosa_juv detections    8685

#  5  2022 cucumaria_frondosa     human_counts 32940
#  6  2022 cucumaria_frondosa     detections   25987
#  7  2022 cucumaria_frondosa_juv human_counts 60714
#  8  2022 cucumaria_frondosa_juv detections    3972

#  9  2023 cucumaria_frondosa     human_counts 43213.
# 10  2023 cucumaria_frondosa     detections   22443
# 11  2023 cucumaria_frondosa_juv human_counts 14640.
# 12  2023 cucumaria_frondosa_juv detections    4531

##########################################################################
# ----------------------------EXPORT DATA --------------------
##########################################################################
head(dplot)
# Need to bind with the original metadata per year


export <- rename(dplot, Station=trial)
head(export)
tail(export)
filter(export, Species=="cucumaria_frondosa") %>% arrange(Station, desc=T)
write.csv(export, "export-2021-2022-2023.csv")

head(dplot)



#########################################################
# Export for Joel (at the ind level) with time stamp
#########################################################
head(d21)
head(d22)
head(d23)

d <- bind_rows(data.frame(select(d21, -min_x), year=2021), data.frame(d22, year=2022), data.frame(mutate(d23, trial=as.character(trial)), year=2023))
head(d)
d <- d %>% select(-dist) %>% rename(ID=Track.id, Station=trial, position.x=max_x)
head(d)
tail(d)

write_csv(d, "SPM_IA_predictions.csv")
#########################################################











# Low counts
ggplot(melt(select(filter(dplot, human_counts < 500), year, human_counts, detections), id.vars=c("year")), aes(x=factor(year), y=value, fill=variable)) + geom_boxplot() + scale_fill_discrete("") + scale_x_discrete("Year") + scale_y_continuous("Counts")

ggplot(filter(dplot, human_counts < 500), aes(human_counts, detections)) + geom_line(aes(group=trial)) + geom_point(aes(color=factor(year))) + geom_smooth(aes(color=factor(year)), method='lm') + theme_light() + xlab("Human count") + ylab("Predicted count") + geom_abline(aes(slope=1, intercept=0), color="black", linetype=3) + scale_color_discrete("YEAR")



# High counts
ggplot(melt(select(filter(dplot, human_counts > 500), year, human_counts, detections), id.vars=c("year")), aes(x=factor(year), y=value, fill=variable)) + geom_boxplot() + scale_fill_discrete("") + scale_x_discrete("Year") + scale_y_continuous("Counts")

ggplot(filter(dplot, human_counts > 500), aes(human_counts, detections)) + geom_line(aes(group=trial)) + geom_point(aes(color=factor(year))) + geom_smooth(aes(color=factor(year)), method='lm') + theme_light() + xlab("Human count") + ylab("Predicted count") + geom_abline(aes(slope=1, intercept=0), color="black", linetype=3) + scale_color_discrete("YEAR")






ggsave("2021 vs 2022 human machine.pdf")



melt(select(d2y, year, human_counts, detections), id.vars=c("year")) %>% group_by(year, variable) %>% summarize(count=sum(value))

melt(select(d2y, year, human_counts, detections), id.vars=c("year")) %>% group_by(year, variable) %>% summarize(count=sum(value)), aes(x=year, y=count, fill=variable)
# 1  2021 human_counts 29189.
# 2  2021 detections   25540
# 3  2022 human_counts 34690
# 4  2022 detections   23356


# 2021 - 12.7% of difference man vs machine
# 2022 - 12.7% of difference man vs machine










export <- left_join(dplot, rename(meta22, trial=station))
filter(export, Species=="cucumaria_frondosa") %>% arrange(trial, desc=T)

write.csv(export, "export.csv")

































































+ geom_errorbar(aes(y=count, xmin=adult_per_station-adult_per_station*min(errorrate$errorrate.min), xmax=adult_per_station+adult_per_station*min(errorrate$errorrate.min)), data=dplotcf) + geom_errorbar(aes(y=count, xmin=adult_per_station-adult_per_station*max(errorrate$errorate.max), xmax=adult_per_station+adult_per_station*max(errorrate$errorate.max)), alpha=0.34, data=dplotcf) + ggtitle(str_c("Th=", th)) + scale_x_continuous(trans='log1p') + scale_y_continuous(trans='log1p') + geom_text(aes(x=adult_per_station, y=count, label=trial, color=Valid), data=dplotcf) + theme_light() + xlab("Human count") + ylab("Predicted count") + geom_abline(aes(slope=1, intercept=0), color="black", linetype=3)
# ggsave("th6.pdf")


ggplot() + geom_point(aes(x=count, y=adult_per_station, color=Valid), data=dplotcf) + theme_light() + xlab("Human count") + ylab("Predicted count") + geom_abline(aes(slope=1, intercept=0), color="black", linetype=3) + geom_errorbar(aes(x=count, ymin=adult_per_station-adult_per_station*min(errorrate$errorrate.min), ymax=adult_per_station+adult_per_station*min(errorrate$errorrate.min)), data=dplotcf) + geom_errorbar(aes(x=count, ymin=adult_per_station-adult_per_station*max(errorrate$errorate.max), ymax=adult_per_station+adult_per_station*max(errorrate$errorate.max)), alpha=0.34, data=dplotcf) + ggtitle(str_c("Th=", th)) +  geom_text_repel(aes(y=adult_per_station, x=count, label=trial, color=Valid), data=dplotcf) + theme_light() + xlab("Human count") + ylab("Predicted count") + geom_abline(aes(slope=1, intercept=0), color="black", linetype=3)





# Filter detections at the bottom of the frames only
# --------------------------------------------------------------------

# Read one track file
stn <- filter(data22, trial==167)
head(stn)

# Check distribution of the counts on the frame
ggplot(stn) + geom_bin2d(aes(x=TL_x, y=TL_y)) + scale_fill_continuous() + scale_y_reverse()

ggplot(filter(stn, (max(TL_y)-min(TL_y)) > dist.min)) + geom_bin2d(aes(x=TL_x, y=TL_y)) + scale_fill_continuous() + scale_y_reverse()



filter(stn, TL_y > 0, Confidence > 0.7) %>% group_by(Species) %>% summarize(count=length(unique(Track.id)))
ggplot(filter(stn, TL_y > 0, Confidence > 0.7, n.frames < 50)) + geom_bin2d(aes(x=TL_x, y=TL_y)) + scale_fill_continuous() + scale_y_reverse()

ggplot(filter(filter(data22, trial%in%unique(data22$trial)[1:16]), TL_y > 0, Confidence > 0.5, n.frames)) + geom_bin2d(aes(x=TL_x, y=TL_y)) + scale_fill_continuous(trans='log1p') + scale_y_reverse() + facet_wrap(~trial)









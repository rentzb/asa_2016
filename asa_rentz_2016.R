### Script Written by Bradley T. Rentz
### November 2016
### 5th Joint Meeting Acoustical Society of America and Acoustical Society of Japan, Honolulu, Hawai'i
### For poster: 2pSC36. The Pohnpeian stop contrast between laminal alveolars and apical dentals involves differences in VOT and F2 locus equation intercepts by Bradley Rentz & Victoria Anderson
### Email: rentzb@hawaii.edu
### rentz.weebly.com

library(dplyr)
library(ggplot2)
library(rstanarm)
library(yarrr)

options (mc.cores=parallel::detectCores ()) # Run on multiple cores

set.seed (3875)

setwd("~/Documents/UH/QP/qp1_data/results") # change wd where files are saved

##############  Locus equations
t_locus <- read.csv("t_d_locus.csv")


# clean data (removed 0 obs) checking for outliers outside of +- 3sd
t_locus_clean <- t_locus %>%
  group_by(speaker,consonant,vowel) %>%
  mutate(sdmeanf2ss = mean(f2ss)+3*sd(f2ss),negsdmeanf2ss=mean(f2ss)-3*sd(f2ss),
         sdmeanf2i=mean(f2i)+3*sd(f2i),negsdmeanf2i=mean(f2i)-3*sd(f2i)) %>%
  filter(f2i<sdmeanf2i)%>%
  filter(f2i>negsdmeanf2i)%>%
  filter(f2ss<sdmeanf2ss)%>%
  filter(f2ss>negsdmeanf2ss)

#subset by consonant t=laminal, d=apical (from PNI orthography)
t_locus_t <- t_locus[t_locus$consonant=="t̻",]
t_locus_d <- t_locus[t_locus$consonant=="t",]

## Bayesian linear hierarchical model (can take a minute or two to run) [laminal]
t_locus.blmer = stan_lmer(f2i ~ f2ss + (1+f2ss|speaker), data=t_locus_t, 
                          prior_intercept = normal(0, 50),
                          prior = normal(0, 2),
                          prior_covariance = decov(regularization = 2), 
                          chains = 4,
                          iter = 2000,adapt_delta=0.9999)

## graph the outputs
draws <- as.data.frame(t_locus.blmer)
colnames(draws)[1:2] <- c("a", "b")

t_locus_plot <- ggplot(t_locus_t, aes(x = f2ss, y = f2i)) + 
  geom_abline(data = draws, aes(intercept = a, slope = b), color = "skyblue", size = 0.2, alpha = 0.1)+
  geom_point(size = 1) + xlim(1000,3000) + ylim(1000,3000) + theme_bw() + xlab("F2 steady state (Hz)")+
  ylab("F2 Initial (Hz)") + ggtitle("Laminal Alveolar Locus Equations") 
t_locus_plot

## Bayesian regression [apical] (takes a minute or two to run)
d_locus.blmer = stan_lmer(f2i ~ f2ss + (1+f2ss|speaker), data=t_locus_d, 
                          prior_intercept = normal(0, 50),
                          prior = normal(0, 2),
                          prior_covariance = decov(regularization = 2), 
                          chains = 4,
                          iter = 2000, adapt_delta = 0.99999 )

# Plot results
draws2 <- as.data.frame(t_locus.blmer)
colnames(draws2)[1:2] <- c("a", "b")

d_locus_plot <- ggplot(t_locus_d, aes(x = f2ss, y = f2i)) + 
  geom_abline(data = draws2, aes(intercept = a, slope = b), color = "seagreen4", size = 0.2, alpha = 0.1)+
  geom_point(size = 1) + xlim(1000,3000) + ylim(1000,3000) + theme_bw() + xlab("F2 steady state (Hz)")+
  ylab("F2 Initial (Hz)") + ggtitle("Apical Dental Locus Equations") 
d_locus_plot


## Pirate plots 
# F2 initial
pirateplot(f2i~vowel+consonant,data=t_locus,inf="hdi",theme=3,hdi.iter=50000,avg.line.fun = median,ylab="F2 Initial (Hz)",gl.col = "white",ylim=c(1200,3200))
# F2 steady state
pirateplot(f2ss~vowel+consonant,data=t_locus,inf="hdi",theme=3,hdi.iter=50000,avg.line.fun = median,ylab="F2 Steady State (Hz)",gl.col = "white",ylim=c(1200,3200))


##############  VOT
#import csv
t_vot <- read.csv("t_d_vot.csv")


# Bayesian regression (takes a minute or two to run) [consonant is contrast coded so consonant output is diffence between the two consonants]
t_vot.clean.blmer = stan_lmer(vot*1000 ~ consonant + (1+consonant|speaker_id) + (1+consonant|word), data=t_vot, 
                              prior_intercept = normal(0, 5),
                              prior = normal(0, 2),
                              prior_covariance = decov(regularization = 2), 
                              chains = 4,
                              iter = 2000,
                              adapt_delta=0.999)

summary(t_vot.clean.blmer)

## Plot pirate plots of VOT (*1000 to convert to ms)
pirateplot(vot*1000~consonant2,data=t_vot,inf="hdi",theme=3,hdi.iter=50000,avg.line.fun = median,ylab="VOT (ms)",gl.col = "white",xlab="Consonant",pal="gray",bean.f.col=c("seagreen4","skyblue"))

##############  Spectral Moments
#import csv
frication <- read.csv("t_d_frication.csv")

# clean data (removes 3sigma outlier) [removed 0 obs.]
frication.clean <- frication %>%
  group_by(speaker,consonant)%>%
  mutate(skewsdmean = mean(skewness) + 3*sd(skewness),negskewsdmean = mean(skewness)-3*sd(skewness),cogsdmean = mean(cog) + 3*sd(cog),negcogsdmean = mean(cog)-3*sd(cog))%>%
  filter(cog<cogsdmean)%>%
  filter(cog>negcogsdmean) %>%
  filter(skewness<skewsdmean)%>%
  filter(skewness>negskewsdmean)

## Bayesian regression for center of gravity
frication.clean.cog.blmer = stan_lmer(cog ~ consonant + (1+consonant|speaker) + (1+consonant|word), data=frication.clean, 
                                      prior_intercept = normal(0, 100),
                                      prior = normal(0, 50),
                                      prior_covariance = decov(regularization = 2), 
                                      chains = 4,
                                      iter = 2000,
                                      adapt_delta=0.999)

plot(frication.clean.cog.blmer,pars=c("(Intercept)","consonantt","consonantt̻"))

# saves regression output to text file in wd
sink("sink-summary-blmer-cog.txt")
print(summary(frication.clean.cog.blmer),digits=4)
sink()

## regression for skewness
frication.clean.skew.blmer = stan_lmer(skewness ~ consonant + (1+consonant|speaker) + (1+consonant|word), data=frication.clean, 
                                       prior_intercept = normal(0, 1),
                                       prior = normal(0, 2),
                                       prior_covariance = decov(regularization = 2), 
                                       chains = 4,
                                       iter = 2000,
                                       adapt_delta=0.999)

## saves regression output as text file in wd
sink("sink-summary-blmer-skew.txt")
print(summary(frication.clean.skew.blmer),digits=3)
sink()

## Pirate plots of data
# skewness
pirateplot(skewness~consonant,data=frication,inf="hdi",theme=3,hdi.iter=50000,avg.line.fun = median,ylab="Skewness",gl.col = "white",xlab="Consonant",pal="gray",bean.f.col=c("purple","seagreen4","skyblue"))
# center of gravity
pirateplot(cog~consonant,data=frication,inf="hdi",theme=3,hdi.iter=50000,avg.line.fun = median,ylab="Center of Gravity (Hz)",gl.col = "white",xlab="Consonant",pal="gray",bean.f.col=c("purple","seagreen4","skyblue"))

#--------------------------------------------------------------
#Ben Neely
#08/08/2024
#Detection probability for Flathead Catfish in small impoundments
#--------------------------------------------------------------

## Clear R
cat("\014")  
rm(list=ls())

## Install and load packages
## Checks if package is installed, installs if not, activates for current session
if("FSA" %in% rownames(installed.packages()) == FALSE) {install.packages("FSA")}
library(FSA)

if("rio" %in% rownames(installed.packages()) == FALSE) {install.packages("rio")}
library(rio)

if("patchwork" %in% rownames(installed.packages()) == FALSE) {install.packages("patchwork")}
library(patchwork)

if("MuMIn" %in% rownames(installed.packages()) == FALSE) {install.packages("MuMIn")}
library(MuMIn)

if("unmarked" %in% rownames(installed.packages()) == FALSE) {install.packages("unmarked")}
library(unmarked)

if("tidyverse" %in% rownames(installed.packages()) == FALSE) {install.packages("tidyverse")}
library(tidyverse)

## Set ggplot theme
pubtheme=theme_classic()+
  theme(panel.grid=element_blank(), 
        panel.background=element_blank(),
        plot.background=element_blank(),
        panel.border=element_rect(color="black",fill="transparent",linewidth=1),
        axis.title=element_text(size=22,color="black",face="bold"),
        axis.text=element_text(size=18,color="black"),
        legend.position=c(0.0001,0.9999),
        legend.justification=c("left","top"),
        legend.direction="horizontal",
        legend.background=element_rect(fill="transparent"),
        legend.title=element_text(size=16),
        legend.text=element_text(size=14),
        strip.text=element_text(size=16))
options(scipen=999)



## Set working directory
setwd("C:/Users/Ben.Neely/OneDrive - State of Kansas, OITS/Documents/Active manuscripts/FCF project/FCF JFWM submission/FCF JFWM third submission")

## Read in data with import
dat=import("Data S1.csv")

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Run everything bracketed here in R Console if it aborts in R Studio

## Occupancy Models
## Create unmarkedFrameOccu object for model fitting
## Set site covariates as prop_samp, mean_depth,med_secchi,pop_dens, and ha
## observation covariates as temp and wind
catfish_occu_dat=unmarkedFrameOccu(y=dat[,3:11],
                                   siteCovs=data.frame(impd=dat[,1],
                                                       prop_samp=dat[,12],
                                                       mean_depth=dat[,13]),
                                   obsCovs=list(temp=dat[,14:22]))

## Look at summary
summary(catfish_occu_dat)

## Scale observation covariate (temp)
catfish_occu_dat@obsCovs$temp=scale(catfish_occu_dat@obsCovs$temp)

## Detection probability without covariates
dp=occu(~1
        ~1,
        catfish_occu_dat)

summary(dp)
backTransform(dp,type="det") ##detection probability

################################################################################
## Fit full model that includes three fixed effects and impoundment as random effect
gc()
mod_full=occu(~prop_samp+mean_depth+temp+(1|impd)
              ~1,
              catfish_occu_dat)


## Use MuMIn to fit all subsets and rank by AICc
modlist=dredge(mod_full,rank="AICc")
(modlist1=modlist%>%
  as_tibble())

## Export model list
#export(modlist,"modlist.csv")

## Extract coefficients with SE for each model
mod1=occu(~prop_samp+mean_depth+temp+(1|impd)
          ~1,
          catfish_occu_dat)
summary(mod1)

mod2=occu(~mean_depth+temp+(1|impd)
          ~1,
          catfish_occu_dat)
summary(mod2)

mod3=occu(~prop_samp+mean_depth+(1|impd)
          ~1,
          catfish_occu_dat)
summary(mod3)

mod4=occu(~prop_samp+temp+(1|impd)
          ~1,
          catfish_occu_dat)
summary(mod4)

mod5=occu(~temp+(1|impd)
          ~1,
          catfish_occu_dat)
summary(mod5)

mod6=occu(~mean_depth+(1|impd)
          ~1,
          catfish_occu_dat)
summary(mod6)

mod7=occu(~prop_samp+(1|impd)
          ~1,
          catfish_occu_dat)
summary(mod7)

mod8=occu(~prop_samp+mean_depth+temp+(1|impd)
          ~1,
          catfish_occu_dat)
summary(mod8)

mod9=occu(~(1|impd)
          ~1,
          catfish_occu_dat)
summary(mod9)

################################################################################
## Organize data to plot model predictions based on most parsimonous model
## Range of mean depth, proportion sampled, and temp
depth_range=range(dat$mean_depth)
propsamp_range=range(dat$prop_samp)
temp_range=range(catfish_occu_dat@obsCovs$temp)

## Create vectors of 10 values split between min and max observations
depths=seq(2,6,1)
propsamps=seq(propsamp_range[1],propsamp_range[2],length.out=10)
temps=seq(temp_range[1],temp_range[2],length.out=10)
impd_unique=unique(dat$impd)

## Create all possible combinations of variable values
tmp=expand.grid(impd=impd_unique,
                mean_depth=depths,
                prop_samp=propsamps,
                temp=temps)

## Estimate detection probability from all combinations of variable values
tmp1=predict(mod_full,type="det",newdata=tmp,appendData=T)

## Get scaled temperature back to original values and create data for plotting
temp_sd=sd(c(dat$temp1,dat$temp2,dat$temp3,dat$temp4,dat$temp5,
             dat$temp6,dat$temp7,dat$temp8,dat$temp9))
temp_mean=mean(c(dat$temp1,dat$temp2,dat$temp3,dat$temp4,dat$temp5,
                  dat$temp6,dat$temp7,dat$temp8,dat$temp9))

## Combine data for plotting
plotdat=tmp1%>%
  mutate(temp_orig=temp_mean+(temp*temp_sd))%>%
  select(det_prob=Predicted,se=SE,lower,upper,mean_depth,prop_samp,temp_orig)

## Export plotting data so I can bring it back into RStudio
export(plotdat,"plotdat.csv")

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Plot model predictions
## Bring in data from above -- Can just call the object if not done in R Console
plotdat_line=import("plotdat.csv")

## Get temperature quantiles for a four panel graph
quantile(c(dat$temp1,dat$temp2,dat$temp3,dat$temp4,dat$temp5,
           dat$temp6,dat$temp7,dat$temp8,dat$temp9))

## Set up facet labels based on quantile values
lab1="24.44 to 26.66\u00B0C"
lab2="26.67 to 27.77\u00B0C"
lab3="27.78 to 29.44\u00B0C"
lab4="29.45 to 31.11\u00B0C"

## Organize data a bit to summarize and look nice for plotting
plotdat_line1=plotdat_line%>%
  mutate(temp_quant=case_when(temp_orig<26.7 ~ lab1,
                              temp_orig>=26.7 & temp_orig<27.8 ~ lab2,
                              temp_orig>=27.8 & temp_orig<29.4 ~ lab3,
                              temp_orig>=29.4 ~ lab4),
         temp_quant=factor(temp_quant,levels=c(lab1,lab2,lab3,lab4)))%>%
  group_by(mean_depth,prop_samp,temp_quant)%>%
  summarize(detprob=mean(det_prob),
            detprob_se=sd(det_prob)/sqrt(n()))%>%
  ungroup()%>%
  mutate(mean_depth=factor(mean_depth))%>%
  filter(mean_depth==2|mean_depth==3|mean_depth==4|mean_depth==5|mean_depth==6)

## Create plot
ggplot(plotdat_line1)+
  geom_line(aes(x=prop_samp,y=detprob,color=mean_depth),linewidth=2)+
  geom_ribbon(aes(x=prop_samp,ymin=detprob-detprob_se,ymax=detprob+detprob_se,fill=mean_depth),alpha=0.5)+
  scale_color_viridis_d(name="Mean depth (m)")+
  scale_fill_viridis_d(name="Mean depth (m)")+
  scale_x_continuous(seq(0.3,0.8,0.1),
                     name="Proportion of impoundment sampled")+
  scale_y_continuous(seq(0.1,0.9,0.1),
                     name="Detection probability")+
  facet_grid(~temp_quant)+
  pubtheme+
  theme(legend.direction="vertical",
        legend.title=element_text(size=20),
        legend.text=element_text(size=18))

## Export plot
ggsave(plot=last_plot(),"Figure 1.png",width=14,height=7,bg="white")
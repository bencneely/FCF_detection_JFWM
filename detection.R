#--------------------------------------------------------------
#Ben Neely
#04/10/2024
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
        panel.border=element_rect(fill="transparent"),
        axis.title=element_text(size=22,color="black",face="bold"),
        axis.text=element_text(size=18,color="black"),
        legend.position="none")
options(scipen=999)

## Set working directory
setwd("C:/Users/Ben.Neely/OneDrive - State of Kansas, OITS/Desktop/FCF JFWM detection/")

## Read in data with import
dat=import("occudat.csv")

################################################################################
##Occupancy Models
## Create unmarkedFrameOccu object for model fitting
## Set site covariate as impoundment and observation covariates as temp, wind, and day of year
catfish_occu_dat=unmarkedFrameOccu(y=dat[,3:11],
                                   siteCovs=data.frame(impd=factor(dat[,1])),
                                   obsCovs=list(temp=dat[,13:21],
                                                wind=dat[,22:30],
                                                doy=dat[,31:39]))

## Look at summary
summary(catfish_occu_dat)

## Scale observation covariates (temp, wind, and day of year)
#obsCovs(catfish_occu_dat)=scale(obsCovs(catfish_occu_dat))

catfish_occu_dat@obsCovs$doy=scale(catfish_occu_dat@obsCovs$doy)
catfish_occu_dat@obsCovs$temp=scale(catfish_occu_dat@obsCovs$temp)
catfish_occu_dat@obsCovs$wind=scale(catfish_occu_dat@obsCovs$wind)

################################################################################
## Fit full model that allows detection to deviate by day, impd, temp, and wind
mod_full=occu(~doy+impd+temp+wind ~1,catfish_occu_dat)

## Use MuMIn to fit all subsets and rank by AICc
modlist=dredge(mod_full,rank="AICc")
(modlist1=modlist%>%
  as_tibble())
export(modlist1,"mods_out.csv")

## Only model within 2 dAICc is the doy + impd + temp model
mod1=occu(~doy+impd+temp ~1,catfish_occu_dat)
summary(mod1)
backTransform(mod1,type="state") ##occupancy probability

## Estimate detection prob for each impoundment allowing impd, doy, and temp to differ
## coefficients are in alphabetical order, in this example...
## crsl, doy, gesl, lvsl, mcla, mgsl, nosl, sgca, wlsl, temp
all_det=backTransform(linearComb(mod1,coefficients=c(0,1,0,0,0,0,0,0,0,1),type="det")) ##all
crsl_det=backTransform(linearComb(mod1,coefficients=c(1,1,0,0,0,0,0,0,0,1),type="det")) ##crsl
gesl_det=backTransform(linearComb(mod1,coefficients=c(0,1,1,0,0,0,0,0,0,1),type="det")) ##gesl
lvsl_det=backTransform(linearComb(mod1,coefficients=c(0,1,0,1,0,0,0,0,0,1),type="det")) ##lvsl
mcla_det=backTransform(linearComb(mod1,coefficients=c(0,1,0,0,1,0,0,0,0,1),type="det")) ##mcla
mgsl_det=backTransform(linearComb(mod1,coefficients=c(0,1,0,0,0,1,0,0,0,1),type="det")) ##mgsl
nosl_det=backTransform(linearComb(mod1,coefficients=c(0,1,0,0,0,0,1,0,0,1),type="det")) ##nosl
sgca_det=backTransform(linearComb(mod1,coefficients=c(0,1,0,0,0,0,0,1,0,1),type="det")) ##sgca
wlsl_det=backTransform(linearComb(mod1,coefficients=c(0,1,0,0,0,0,0,0,1,1),type="det")) ##wlsl

impds=c("all","crsl","gesl","lvsl","mcla","mgsl","nosl","sgca","wlsl")
det=c(all_det@estimate,crsl_det@estimate,gesl_det@estimate,lvsl_det@estimate,mcla_det@estimate,
      mgsl_det@estimate,nosl_det@estimate,sgca_det@estimate,wlsl_det@estimate)
stderr=c(SE(all_det),SE(crsl_det),SE(gesl_det),SE(lvsl_det),SE(mcla_det),
         SE(mgsl_det),SE(nosl_det),SE(sgca_det),SE(wlsl_det))

detprobs=bind_cols(impd=impds,det_prob=det,det_se=stderr)

################################################################################
## Create plots of detection probability for each pop
## Set min, mean, and max doy and temp
min_doy=round(min(catfish_occu_dat@obsCovs$doy),1)
med_doy=round(median(catfish_occu_dat@obsCovs$doy),1)
max_doy=round(max(catfish_occu_dat@obsCovs$doy),1)

min_temp=round(min(catfish_occu_dat@obsCovs$temp),1)
max_temp=round(max(catfish_occu_dat@obsCovs$temp),1)

###########################################################  
## sgca
## Set up new data frame with fixed doy and temp to estimate detection prob at each
sgca=expand.grid(impd="sgca",
                 doy=seq(from=min_doy,to=max_doy,by=0.1),
                 temp = seq(from=min_temp,to=max_temp,by=0.01))

## Predict detection prob at each measure of doy and temp (scaled)                 
sgca_pred=predict(mod1,type ="det",newdata=sgca,appendData = TRUE)

## Plot detection probability at sgca for earliest, mean, and latest doy at temp
sgca_plot=ggplot()+
  geom_line(subset(sgca_pred,doy==min_doy),mapping=aes(x=temp,y=Predicted),color="#440154",linewidth=2)+
  geom_line(subset(sgca_pred,doy==med_doy),mapping=aes(x=temp,y=Predicted),color="#238A8D",linewidth=2)+
  geom_line(subset(sgca_pred,doy==max_doy),mapping=aes(x=temp,y=Predicted),color="#FDE725",linewidth=2)+
  geom_ribbon(subset(sgca_pred,doy==min_doy),mapping=aes(x=temp,ymin=lower,ymax=upper),fill="#440154",alpha=0.2)+
  geom_ribbon(subset(sgca_pred,doy==med_doy),mapping=aes(x=temp,ymin=lower,ymax=upper),fill="#238A8D",alpha=0.2)+
  geom_ribbon(subset(sgca_pred,doy==max_doy),mapping=aes(x=temp,ymin=lower,ymax=upper),fill="#FDE725",alpha=0.2)+
  scale_x_continuous(limits=c(min_temp-0.01,max_temp+0.02),
                     breaks=seq(-2,2,0.5),
                     name="")+
  scale_y_continuous(limits=c(0,1),
                     breaks=seq(0,1,0.1),
                     name="Detection probability")+
  annotate("text",x=min_temp,y=1,label="Afton",hjust=0,vjust=1,size=8)+
  coord_cartesian(xlim=c(min_temp-0.01,max_temp+0.1),
                  ylim=c(-0.01,1.01),
                  expand=F)+
  pubtheme

###########################################################  
## crsl
## Set up new data frame with fixed doy and temp to estimate detection prob at each
crsl=expand.grid(impd="crsl",
                 doy=seq(from=min_doy,to=max_doy,by=0.1),
                 temp = seq(from=min_temp,to=max_temp,by=0.01))

## Predict detection prob at each measure of doy and temp (scaled)                 
crsl_pred=predict(mod1,type ="det",newdata=crsl,appendData = TRUE)

## Plot detection probability at crsl for earliest, mean, and latest doy at temp
crsl_plot=ggplot()+
  geom_line(subset(crsl_pred,doy==min_doy),mapping=aes(x=temp,y=Predicted),color="#440154",linewidth=2)+
  geom_line(subset(crsl_pred,doy==med_doy),mapping=aes(x=temp,y=Predicted),color="#238A8D",linewidth=2)+
  geom_line(subset(crsl_pred,doy==max_doy),mapping=aes(x=temp,y=Predicted),color="#FDE725",linewidth=2)+
  geom_ribbon(subset(crsl_pred,doy==min_doy),mapping=aes(x=temp,ymin=lower,ymax=upper),fill="#440154",alpha=0.2)+
  geom_ribbon(subset(crsl_pred,doy==med_doy),mapping=aes(x=temp,ymin=lower,ymax=upper),fill="#238A8D",alpha=0.2)+
  geom_ribbon(subset(crsl_pred,doy==max_doy),mapping=aes(x=temp,ymin=lower,ymax=upper),fill="#FDE725",alpha=0.2)+
  scale_x_continuous(limits=c(min_temp-0.01,max_temp+0.02),
                     breaks=seq(-2,2,0.5),
                     name="")+
  scale_y_continuous(limits=c(0,1),
                     breaks=seq(0,1,0.1),
                     name="Detection probability")+
  annotate("text",x=min_temp,y=1,label="Crawford",hjust=0,vjust=1,size=8)+
  coord_cartesian(xlim=c(min_temp-0.01,max_temp+0.1),
                  ylim=c(-0.01,1.01),
                  expand=F)+
  pubtheme

###########################################################  
## gesl
## Set up new data frame with fixed doy and temp to estimate detection prob at each
gesl=expand.grid(impd="gesl",
                 doy=seq(from=min_doy,to=max_doy,by=0.1),
                 temp = seq(from=min_temp,to=max_temp,by=0.01))

## Predict detection prob at each measure of doy and temp (scaled)                 
gesl_pred=predict(mod1,type ="det",newdata=gesl,appendData = TRUE)

## Plot detection probability at gesl for earliest, mean, and latest doy at temp
gesl_plot=ggplot()+
  geom_line(subset(gesl_pred,doy==min_doy),mapping=aes(x=temp,y=Predicted),color="#440154",linewidth=2)+
  geom_line(subset(gesl_pred,doy==med_doy),mapping=aes(x=temp,y=Predicted),color="#238A8D",linewidth=2)+
  geom_line(subset(gesl_pred,doy==max_doy),mapping=aes(x=temp,y=Predicted),color="#FDE725",linewidth=2)+
  geom_ribbon(subset(gesl_pred,doy==min_doy),mapping=aes(x=temp,ymin=lower,ymax=upper),fill="#440154",alpha=0.2)+
  geom_ribbon(subset(gesl_pred,doy==med_doy),mapping=aes(x=temp,ymin=lower,ymax=upper),fill="#238A8D",alpha=0.2)+
  geom_ribbon(subset(gesl_pred,doy==max_doy),mapping=aes(x=temp,ymin=lower,ymax=upper),fill="#FDE725",alpha=0.2)+
  scale_x_continuous(limits=c(min_temp-0.01,max_temp+0.02),
                     breaks=seq(-2,2,0.5),
                     name="")+
  scale_y_continuous(limits=c(0,1),
                     breaks=seq(0,1,0.1),
                     name="Detection probability")+
  annotate("text",x=min_temp,y=1,label="Geary",hjust=0,vjust=1,size=8)+
  coord_cartesian(xlim=c(min_temp-0.01,max_temp+0.1),
                  ylim=c(-0.01,1.01),
                  expand=F)+
  pubtheme

###########################################################  
## lvsl
## Set up new data frame with fixed doy and temp to estimate detection prob at each
lvsl=expand.grid(impd="lvsl",
                 doy=seq(from=min_doy,to=max_doy,by=0.1),
                 temp = seq(from=min_temp,to=max_temp,by=0.01))

## Predict detection prob at each measure of doy and temp (scaled)                 
lvsl_pred=predict(mod1,type ="det",newdata=lvsl,appendData = TRUE)

## Plot detection probability at lvsl for earliest, mean, and latest doy at temp
lvsl_plot=ggplot()+
  geom_line(subset(lvsl_pred,doy==min_doy),mapping=aes(x=temp,y=Predicted),color="#440154",linewidth=2)+
  geom_line(subset(lvsl_pred,doy==med_doy),mapping=aes(x=temp,y=Predicted),color="#238A8D",linewidth=2)+
  geom_line(subset(lvsl_pred,doy==max_doy),mapping=aes(x=temp,y=Predicted),color="#FDE725",linewidth=2)+
  geom_ribbon(subset(lvsl_pred,doy==min_doy),mapping=aes(x=temp,ymin=lower,ymax=upper),fill="#440154",alpha=0.2)+
  geom_ribbon(subset(lvsl_pred,doy==med_doy),mapping=aes(x=temp,ymin=lower,ymax=upper),fill="#238A8D",alpha=0.2)+
  geom_ribbon(subset(lvsl_pred,doy==max_doy),mapping=aes(x=temp,ymin=lower,ymax=upper),fill="#FDE725",alpha=0.2)+
  scale_x_continuous(limits=c(min_temp-0.01,max_temp+0.02),
                     breaks=seq(-2,2,0.5),
                     name="Scaled temperature")+
  scale_y_continuous(limits=c(0,1),
                     breaks=seq(0,1,0.1),
                     name="Detection probability")+
  annotate("text",x=min_temp,y=1,label="Leavenworth",hjust=0,vjust=1,size=8)+
  coord_cartesian(xlim=c(min_temp-0.01,max_temp+0.1),
                  ylim=c(-0.01,1.01),
                  expand=F)+
  pubtheme

###########################################################  
## mcla
## Set up new data frame with fixed doy and temp to estimate detection prob at each
mcla=expand.grid(impd="mcla",
                 doy=seq(from=min_doy,to=max_doy,by=0.1),
                 temp = seq(from=min_temp,to=max_temp,by=0.01))

## Predict detection prob at each measure of doy and temp (scaled)                 
mcla_pred=predict(mod1,type ="det",newdata=mcla,appendData = TRUE)

## Plot detection probability at mcla for earliest, mean, and latest doy at temp
mcla_plot=ggplot()+
  geom_line(subset(mcla_pred,doy==min_doy),mapping=aes(x=temp,y=Predicted),color="#440154",linewidth=2)+
  geom_line(subset(mcla_pred,doy==med_doy),mapping=aes(x=temp,y=Predicted),color="#238A8D",linewidth=2)+
  geom_line(subset(mcla_pred,doy==max_doy),mapping=aes(x=temp,y=Predicted),color="#FDE725",linewidth=2)+
  geom_ribbon(subset(mcla_pred,doy==min_doy),mapping=aes(x=temp,ymin=lower,ymax=upper),fill="#440154",alpha=0.2)+
  geom_ribbon(subset(mcla_pred,doy==med_doy),mapping=aes(x=temp,ymin=lower,ymax=upper),fill="#238A8D",alpha=0.2)+
  geom_ribbon(subset(mcla_pred,doy==max_doy),mapping=aes(x=temp,ymin=lower,ymax=upper),fill="#FDE725",alpha=0.2)+
  scale_x_continuous(limits=c(min_temp-0.01,max_temp+0.02),
                     breaks=seq(-2,2,0.5),
                     name="")+
  scale_y_continuous(limits=c(0,1),
                     breaks=seq(0,1,0.1),
                     name="")+
  annotate("text",x=min_temp,y=1,label="Middle Creek",hjust=0,vjust=1,size=8)+
  coord_cartesian(xlim=c(min_temp-0.01,max_temp+0.1),
                  ylim=c(-0.01,1.01),
                  expand=F)+
  pubtheme

###########################################################  
## mgsl
## Set up new data frame with fixed doy and temp to estimate detection prob at each
mgsl=expand.grid(impd="mgsl",
                 doy=seq(from=min_doy,to=max_doy,by=0.1),
                 temp = seq(from=min_temp,to=max_temp,by=0.01))

## Predict detection prob at each measure of doy and temp (scaled)                 
mgsl_pred=predict(mod1,type ="det",newdata=mgsl,appendData = TRUE)

## Plot detection probability at mgsl for earliest, mean, and latest doy at temp
mgsl_plot=ggplot()+
  geom_line(subset(mgsl_pred,doy==min_doy),mapping=aes(x=temp,y=Predicted),color="#440154",linewidth=2)+
  geom_line(subset(mgsl_pred,doy==med_doy),mapping=aes(x=temp,y=Predicted),color="#238A8D",linewidth=2)+
  geom_line(subset(mgsl_pred,doy==max_doy),mapping=aes(x=temp,y=Predicted),color="#FDE725",linewidth=2)+
  geom_ribbon(subset(mgsl_pred,doy==min_doy),mapping=aes(x=temp,ymin=lower,ymax=upper),fill="#440154",alpha=0.2)+
  geom_ribbon(subset(mgsl_pred,doy==med_doy),mapping=aes(x=temp,ymin=lower,ymax=upper),fill="#238A8D",alpha=0.2)+
  geom_ribbon(subset(mgsl_pred,doy==max_doy),mapping=aes(x=temp,ymin=lower,ymax=upper),fill="#FDE725",alpha=0.2)+
  scale_x_continuous(limits=c(min_temp-0.01,max_temp+0.02),
                     breaks=seq(-2,2,0.5),
                     name="")+
  scale_y_continuous(limits=c(0,1),
                     breaks=seq(0,1,0.1),
                     name="")+
  annotate("text",x=min_temp,y=1,label="Montgomery",hjust=0,vjust=1,size=8)+
  coord_cartesian(xlim=c(min_temp-0.01,max_temp+0.1),
                  ylim=c(-0.01,1.01),
                  expand=F)+
  pubtheme

###########################################################  
## nosl
## Set up new data frame with fixed doy and temp to estimate detection prob at each
nosl=expand.grid(impd="nosl",
                 doy=seq(from=min_doy,to=max_doy,by=0.1),
                 temp = seq(from=min_temp,to=max_temp,by=0.01))

## Predict detection prob at each measure of doy and temp (scaled)                 
nosl_pred=predict(mod1,type ="det",newdata=nosl,appendData = TRUE)

## Plot detection probability at nosl for earliest, mean, and latest doy at temp
nosl_plot=ggplot()+
  geom_line(subset(nosl_pred,doy==min_doy),mapping=aes(x=temp,y=Predicted),color="#440154",linewidth=2)+
  geom_line(subset(nosl_pred,doy==med_doy),mapping=aes(x=temp,y=Predicted),color="#238A8D",linewidth=2)+
  geom_line(subset(nosl_pred,doy==max_doy),mapping=aes(x=temp,y=Predicted),color="#FDE725",linewidth=2)+
  geom_ribbon(subset(nosl_pred,doy==min_doy),mapping=aes(x=temp,ymin=lower,ymax=upper),fill="#440154",alpha=0.2)+
  geom_ribbon(subset(nosl_pred,doy==med_doy),mapping=aes(x=temp,ymin=lower,ymax=upper),fill="#238A8D",alpha=0.2)+
  geom_ribbon(subset(nosl_pred,doy==max_doy),mapping=aes(x=temp,ymin=lower,ymax=upper),fill="#FDE725",alpha=0.2)+
  scale_x_continuous(limits=c(min_temp-0.01,max_temp+0.02),
                     breaks=seq(-2,2,0.5),
                     name="")+
  scale_y_continuous(limits=c(0,1),
                     breaks=seq(0,1,0.1),
                     name="")+
  annotate("text",x=min_temp,y=1,label="Neosho",hjust=0,vjust=1,size=8)+
  coord_cartesian(xlim=c(min_temp-0.01,max_temp+0.1),
                  ylim=c(-0.01,1.01),
                  expand=F)+
  pubtheme

###########################################################  
## wlsl
## Set up new data frame with fixed doy and temp to estimate detection prob at each
wlsl=expand.grid(impd="wlsl",
                 doy=seq(from=min_doy,to=max_doy,by=0.1),
                 temp = seq(from=min_temp,to=max_temp,by=0.01))

## Predict detection prob at each measure of doy and temp (scaled)                 
wlsl_pred=predict(mod1,type ="det",newdata=wlsl,appendData = TRUE)

## Create sham data for creating a legend
wlsl_sham=bind_cols(Predicted=c(2,3,4),
                    lower=c(2,3,4),
                    upper=c(2,3,4),
                    temp=c(2,3,4),
                    a=factor(c("First sample date",
                               "Median sample date",
                               "Last sample date"),
                             levels=c("First sample date",
                                      "Median sample date",
                                      "Last sample date")))

## Plot detection probability at wlsl for earliest, mean, and latest doy at temp
wlsl_plot=ggplot()+
  geom_line(subset(wlsl_pred,doy==min_doy),mapping=aes(x=temp,y=Predicted),color="#440154",linewidth=2)+
  geom_line(subset(wlsl_pred,doy==med_doy),mapping=aes(x=temp,y=Predicted),color="#238A8D",linewidth=2)+
  geom_line(subset(wlsl_pred,doy==max_doy),mapping=aes(x=temp,y=Predicted),color="#FDE725",linewidth=2)+
  geom_ribbon(subset(wlsl_pred,doy==min_doy),mapping=aes(x=temp,ymin=lower,ymax=upper),fill="#440154",alpha=0.2)+
  geom_ribbon(subset(wlsl_pred,doy==med_doy),mapping=aes(x=temp,ymin=lower,ymax=upper),fill="#238A8D",alpha=0.2)+
  geom_ribbon(subset(wlsl_pred,doy==max_doy),mapping=aes(x=temp,ymin=lower,ymax=upper),fill="#FDE725",alpha=0.2)+
  
  geom_line(wlsl_sham,mapping=aes(x=temp,y=Predicted,color=a),linewidth=2)+
  geom_ribbon(wlsl_sham,mapping=aes(x=temp,ymin=lower,ymax=upper,fill=a),alpha=0.2)+
  scale_color_manual(values=c("#440154","#238A8D","#FDE725"))+
  scale_fill_manual(values=c("#440154","#238A8D","#FDE725"))+
  
  scale_x_continuous(limits=c(min_temp-0.01,max_temp+0.02),
                     breaks=seq(-2,2,0.5),
                     name="Scaled temperature")+
  scale_y_continuous(limits=c(0,1),
                     breaks=seq(0,1,0.1),
                     name="")+
  annotate("text",x=min_temp,y=1,label="Wilson",hjust=0,vjust=1,size=8)+
  coord_cartesian(xlim=c(min_temp-0.01,max_temp+0.1),
                  ylim=c(-0.01,1.01),
                  expand=F)+
  pubtheme+
  theme(legend.position=c(1,0.01),
        legend.justification=c("right","bottom"),
        legend.title=element_blank(),
        legend.text=element_text(size=16),
        legend.background=element_blank(),
        legend.key.width=unit(2,"cm"))

################################################################################
## Create plot
detplots=(sgca_plot/crsl_plot/gesl_plot/lvsl_plot)|(mcla_plot/mgsl_plot/nosl_plot/wlsl_plot)
ggsave(plot=detplots,"detprobs.png",width=12,height=16,units="in",bg="white")

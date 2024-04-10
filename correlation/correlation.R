#--------------------------------------------------------------
#Ben Neely
#04/10/2024
#Look for correlation in survey covariates
#--------------------------------------------------------------

## Clear R
cat("\014")  
rm(list=ls())

## Install and load packages
## Checks if package is installed, installs if not, activates for current session
if("rio" %in% rownames(installed.packages()) == FALSE) {install.packages("rio")}
library(rio)

if("tidyverse" %in% rownames(installed.packages()) == FALSE) {install.packages("tidyverse")}
library(tidyverse)

## Set working directory
setwd("C:/Users/Ben.Neely/OneDrive - State of Kansas, OITS/Desktop/FCF JFWM detection/")

## Read in data with import
dat=import("correlation/corrdat.csv")

############################################################
## Examine data
## Temp vs Wind
ggplot(dat,aes(x=temp,y=wind,color=impd))+
  geom_point()+
  geom_smooth(method="lm",se=F)

## Temp vs DOY
ggplot(dat,aes(x=temp,y=doy,color=impd))+
  geom_point()+
  geom_smooth(method="lm",se=F)

## Wind vs DOY
ggplot(dat,aes(x=wind,y=doy,color=impd))+
  geom_point()+
  geom_smooth(method="lm",se=F)

############################################################
## Look at correlation coefficients
## Temp vs Wind
cor(dat$temp,dat$wind,method="pearson")
#r=-0.086

## Temp vs DOY
cor(dat$temp,dat$doy,method="pearson")
#r=0.134

## Wind vs DOY
cor(dat$wind,dat$doy,method="pearson")
#r=0.285
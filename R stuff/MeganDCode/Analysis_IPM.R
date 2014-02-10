#GET AND ORGANIZE THE DATA

rm(list=ls(all=TRUE)) ## Clear all

setwd("~/Grad school/CAMEO/R stuff") #set working directory to location of this script & data file

getwd() #check working directory is set properly

#read in data
data<-read.table("sortdata_10Nov.txt", sep="\t", header=TRUE)#read in data

#create site names
sitenametable<-c("HonN",
                 "HonS",
                 "Molo",
                 "Puam",
                 "Kahe",
                 "KanB",
                 "KanP",
                 "Maal",
                 "Olow",
                 "Papa")
sitetable<-c("Honolua",
             "Honolua",
             "Molokini",
             "Puamana",
             "Kahekili",
             "Kanahena",
             "Kanahena",
             "Maalaea",
             "Olowalu",
             "Papaula")
data$site<-sitetable[charmatch(data$sitename,sitenametable)]
data$site<-factor(data$site)
data$quad<-factor(paste("Q",data$quad,sep=""))

#create transect names
data$trans <- as.character(data$depthname)
data <- within(data,trans[sitename=="HonS"]<-"S")
data <- within(data,trans[sitename=="HonN"]<-"N")
data$trans <-factor(data$trans)

#create unique names for transect, quad, and colony
data$utrans<-factor(paste(data$site,data$trans,sep="."))
data$uquad <-factor(paste(data$utrans,data$quad,sep="."))
data$ucolony <-factor(paste(data$utrans,data$quad,data$species,data$colony,sep="."))

#make year a factor
data$fyear_t<-factor(data$year_t)
data$fyear_t1<-factor(data$year_t1)                          

#create some additional variables using logical conditions on fate
data$surv = NA   #SURVIVAL
  data <- within(data,surv[fate_t1!="D"]<-1) #if you're not dead, you survived'
  data <- within(data,surv[fate_t1=="D"]<-0) #if you're dead, you're dead
  data <- within(data,surv[fate_t1=="R"|fate_t1=="N"]<-NA) #don't count if new recruit
data$fiss = NA   #FISSION
  data <- within(data,fiss[fate_t1=="I"|fate_t1=="M"]<-1)
  data <- within(data,fiss[fate_t1!="I" & fate_t1!="M"]<-0)
  data <- within(data,fiss[fate_t1=="R"|fate_t1=="N"]<-NA) #don't count if new recruit
data$fuse = NA   #FUSION
  data <- within(data,fuse[fate_t1=="U"|fate_t1=="M"]<-1)
  data <- within(data,fuse[fate_t1!="U" & fate_t1!="M"]<-0)
  data <- within(data,fuse[fate_t1=="R"|fate_t1=="N"]<-NA) #don't count if new recruit

#create some growth measures
data$lnarea_t<-log(data$area_t)
data$lnarea_t1<-log(data$area_t1)
data$lnperi_t <-log(data$peri_t)
data$lnperi_t1 <-log(data$peri_t1)

data$PAratio<-data$peri_t/data$area_t
data$grArea<-log(data$area_t1/data$area_t)
data$grPeri<-log(data$peri_t1/data$peri_t)
data$grPAR<-log((data$area_t1/data$peri_t1)/(data$area_t/data$peri_t))

#check format of data
  str(data) #check structure of data and variables

#create separate dataframes for three most abundant species
#exclude 2-year observations
MC<-subset(data,(species=="MC" & !skipyear))
PL<-subset(data,(species=="PL" & !skipyear))
PC<-subset(data,(species=="PC" & !skipyear))
#data for growth model, "S" only (excludes fission, fusion)
MCg<-subset(MC,(fate_t1=="S" & utrans!="Puamana.13m"))
PLg<-subset(PL,(fate_t1=="S"))
PCg<-subset(PC,(fate_t1=="S"))

#Model growth with linear mixed-effects model in nlme
library(nlme) 

#Plots of fullest model by utrans
grMC_FF.list<-lmList(grArea ~ lnarea_t*spcover_t*allcover_t*depth|utrans,MCg,na.action=na.pass)
plot(intervals(grMC_FF.list))

#Fullest possible fixed effect model includes:
#indiv area, consp cover, total coral cover, depth, exposure

#model for growth, fullest fixed model 
    #use gls to get REML fit for model comparison
grMC_FF <-  gls(grArea~lnarea_t*spcover_t*allcover_t*depth,method="REML",data=MCg)
grPL_FF <-  gls(grArea~lnarea_t*spcover_t*allcover_t*depth, method="REML",data=PLg)
grPC_FF <-  gls(grArea~lnarea_t*spcover_t*allcover_t*depth, method="REML",data=PCg)
    #use lm to fit and get nice residual plots
grMC_FFlm <-  lm(grArea~lnarea_t*spcover_t*allcover_t*depth,data=MCg)
  E_MC<-rstandard(grMC_FFlm)
  boxplot(E_MC~utrans,data=MCg)
grPL_FFlm <-  lm(grArea~lnarea_t*spcover_t*allcover_t*depth,data=PLg)
  E_PL<-rstandard(grPL_FFlm)
  boxplot(E_PL~utrans,data=PLg)
grPC_FFlm <-  lm(grArea~lnarea_t*spcover_t*allcover_t*depth,data=PCg)
  E_PC<-rstandard(grPC_FFlm)
  boxplot(E_PC~utrans,data=PCg)

  #check residuals for fullest fixed model
  plot(grMC_FFlm,which=1:6,id.n=5,labels.id=paste(MCg$year_t,MCg$ucolony,sep="."))
  plot(grPL_FFlm,which=1:6,id.n=5,labels.id=paste(PLg$year_t,PLg$ucolony,sep="."))
  plot(grPC_FFlm,which=1:6,id.n=5,labels.id=paste(PCg$year_t,PCg$ucolony,sep="."))
  #outliers added to list for checking, but nothing bad

#fit best random effects model
lmc<-lmeControl(niter=10000,msMaxIter=10000)
grMC_RE1 <-lme(grArea~lnarea_t*spcover_t*allcover_t*depth,
               random=~1+spcover_t+allcover_t +lnarea_t|utrans,data=MCg)
grMC_RE2 <-lme(grArea~lnarea_t*spcover_t*allcover_t*depth,
               random=~1+allcover_t +lnarea_t|utrans,data=MCg)
grMC_RE3 <-lme(grArea~lnarea_t*spcover_t*allcover_t*depth,
               random=~1+spcover_t+lnarea_t|utrans,data=MCg, control=lmc)
grMC_RE4 <-lme(grArea~lnarea_t*spcover_t*allcover_t*depth,
               random=~1+spcover_t+allcover_t|utrans,data=MCg)
grMC_RE5 <-lme(grArea~lnarea_t*spcover_t*allcover_t*depth,
               random=~1+spcover_t|utrans,data=MCg)
grMC_RE6 <-lme(grArea~lnarea_t*spcover_t*allcover_t*depth,
               random=~1+allcover_t|utrans,data=MCg)
grMC_RE7 <-lme(grArea~lnarea_t*spcover_t*allcover_t*depth,
               random=~1+lnarea_t|utrans,data=MCg,control=lmc)
grMC_RE8 <-lme(grArea~lnarea_t*spcover_t*allcover_t*depth,
               random=~1|utrans,data=MCg, control=lmc)

anova(grMC_FF,grMC_RE1,grMC_RE2,grMC_RE3,grMC_RE4,grMC_RE5,grMC_RE6,grMC_RE7,grMC_RE8)
anova(grMC_RE5,grMC_RE8)
anova(grMC_RE6,grMC_RE8)
anova(grMC_FF,grMC_RE8)
anova(grMC_FF,grMC_RE6)
#USE 6:  lowest BIC is 8, lowest AIC is 6;6 better than 8 by LRT
#including comparison to FF (no random effects) with LRT, RE6 slightly better p=0.047


grPL_RE1 <-lme(grArea~lnarea_t*spcover_t*allcover_t*depth,
               random=~1+spcover_t+allcover_t +lnarea_t|utrans,data=PLg,control=lmc)
grPL_RE2 <-lme(grArea~lnarea_t*spcover_t*allcover_t*depth,
               random=~1+allcover_t +lnarea_t|utrans,data=PLg,control=lmc)
grPL_RE3 <-lme(grArea~lnarea_t*spcover_t*allcover_t*depth,
               random=~1+spcover_t+lnarea_t|utrans,data=PLg)
grPL_RE4 <-lme(grArea~lnarea_t*spcover_t*allcover_t*depth,
               random=~1+spcover_t+allcover_t|utrans,data=PLg)
grPL_RE5 <-lme(grArea~lnarea_t*spcover_t*allcover_t*depth,
               random=~1+spcover_t|utrans,data=PLg)
grPL_RE6 <-lme(grArea~lnarea_t*spcover_t*allcover_t*depth,
               random=~1+allcover_t|utrans,data=PLg)
grPL_RE7 <-lme(grArea~lnarea_t*spcover_t*allcover_t*depth,
               random=~1+lnarea_t|utrans,data=PLg)
grPL_RE8 <-lme(grArea~lnarea_t*spcover_t*allcover_t*depth,
               random=~1|utrans,data=PLg)
anova(grPL_FF,grPL_RE1,grPL_RE2,grPL_RE3,grPL_RE4,grPL_RE5,grPL_RE6,grPL_RE7,grPL_RE8)
anova(grPL_RE2,grPL_RE8)           
anova(grPL_FF,grPL_RE8)
#USE 8:  2 and 8 have siilar AIC, 8 is best BIC; LRT: 2 not better than 8, p=0.07; 
               #RE8 better than FF by LRT (p=0.043)

grPC_RE1 <-lme(grArea~lnarea_t*spcover_t*allcover_t*depth,
               random=~1+spcover_t+allcover_t +lnarea_t|utrans,data=PCg)
grPC_RE2 <-lme(grArea~lnarea_t*spcover_t*allcover_t*depth,
               random=~1+allcover_t +lnarea_t|utrans,data=PCg)
grPC_RE3 <-lme(grArea~lnarea_t*spcover_t*allcover_t*depth,
               random=~1+spcover_t+lnarea_t|utrans,data=PCg)
grPC_RE4 <-lme(grArea~lnarea_t*spcover_t*allcover_t*depth,
               random=~1+spcover_t+allcover_t|utrans,data=PCg)
grPC_RE5 <-lme(grArea~lnarea_t*spcover_t*allcover_t*depth,
               random=~1+spcover_t|utrans,data=PCg)
grPC_RE6 <-lme(grArea~lnarea_t*spcover_t*allcover_t*depth,
               random=~1+allcover_t|utrans,data=PCg)
grPC_RE7 <-lme(grArea~lnarea_t*spcover_t*allcover_t*depth,
               random=~1+lnarea_t|utrans,data=PCg)
grPC_RE8 <-lme(grArea~lnarea_t*spcover_t*allcover_t*depth,
               random=~1|utrans,data=PCg)
anova(grPC_FF,grPC_RE1,grPC_RE2,grPC_RE3,grPC_RE4,grPC_RE5,grPC_RE6,grPC_RE7,grPC_RE8)
 


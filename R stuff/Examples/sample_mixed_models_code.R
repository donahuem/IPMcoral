rm(list=ls(all=TRUE)) ## Clear all

setwd("~/2011 Fall/IPM group") #set working directory to location of this script & data file
  #In RStudio use Tools>Set working directory>To Source file location
getwd() #check working directory is set properly

#read in data
data<-read.table("sample.data.txt", sep="\t", header=TRUE)#read in data
  #pretend this data comes from three sites, each with individuals in high and low light treatments
  #for part of the examples, i'll pretend that this has a split plot design, with on low-light and one high-light plot per site
#check format of data
    str(data) #check structure of data and variables
    data$tag<-as.factor(data$tag)#if you want to change tag number to factor because the number is not meaningfull
    data$light<-factor(data$light,levels(data$light)[c(2,1)])#factor levels default to alphabetic order, this changes order to low, high
    str(data)


#Modeling growth with linear mixed-effects model in lme4
library(lme4) #for lmer/glmer functions (has better glmer function than nlme package)
mod1<-lmer(size.end~light+size.start+(1|site/light), data=data, REML=T) #intercept can vary with site and with split plot (light)
    mod1 #see summary of results
    ranef(mod1) #eBLUPs
    fixef(mod1)#fixed parameter estimates
    resid(mod1)#residuals
    fitted(mod1)#fits
  #model checks
    qqnorm(resid(mod1))#quantile-quantile plot of residuals
    hist(resid(mod1)) #histogram of residuals
    library(lattice) #for smoother
    xyplot(resid(mod1)~fitted(mod1), #residuals vs fitted with smoother
                panel=function(x,y){panel.xyplot(x,y); panel.loess(x,y, span=0.75)})
    plot(resid(mod1)~mod1@frame$site) #residuals by site
    plot(resid(mod1)~mod1@frame$light) #residuals by light trt
    xyplot(resid(mod1)~mod1@frame$size.start, #residuals by starting size
                panel=function(x,y){panel.xyplot(x,y); panel.loess(x,y, span=0.75)})

#test significance of split plot (light) as a random effect [though may want to leave it is as part of experimental design]
    mod2<-lmer(size.end~light+size.start+(1|site), data=data, REML=T)
      anova(mod1, mod2) #not significant

#test significance of light as fixed effect by LRT (switch to ML for this)
    mod1ML<-lmer(size.end~light+size.start+(1|site/light), data=data, REML=F) #change "REML" to false
    mod3ML<-lmer(size.end~size.start+(1|site/light), data=data, REML=F)
      anova(mod1ML, mod3ML) #p=0.01241

#Modeling survival with generalized linear mixed effects model
s1<-lmer(survived~size.start+light+(1|site/light), data=data, family="binomial")
s1
    xyplot(resid(s1)~fitted(s1), #not so useful for binomial
                panel=function(x,y){panel.xyplot(x,y); panel.loess(x,y, span=0.75)})

#Alternate option: modeling growth in nlme
#Allows specification of variance of residuals and correlation of residuals
detach("package:lme4")#can't have both lme4 and nlme at once
library(nlme)
#model without fixed effects (same as linear regression, but in a form R can compare with mixed effects models)
g1<-gls(size.end~size.start+light, na.action=na.omit, data=data, method="REML")
  plot(resid(g1)~fitted(g1))
#generalized least squares, with variance of residuals proportional to starting size (also could make power or exponential function)
g2<-gls(size.end~size.start+light, na.action=na.omit, weight=varFixed(~size.start), data=data, method="REML")
  plot(resid(g2)~fitted(g2))
  anova(g1, g2) #prefer g1 (makes sense -- no evidence of change in variance with size in this data)
#mixed-effects model with size.start and light as fixed, and site and subplot (light) nested within site as random:
g3<-lme(size.end~size.start+light, random=~1|site/light, na.action=na.omit, data=data, method="REML")
    ranef(g3)
    fixef(g3) #estimates essentially the same as mod1, slightly dif from g1
    plot(resid(g3)~data$site[complete.cases(data$size.end)]) #complete.cases removes NAs
    anova(g1, g3) #not significant variation by site/split plot
#Two ways of specifying same model, with random intercept for site:
    g4<-gls(size.end~size.start+light, correlation=corCompSymm(form=~1|site), na.action=na.omit, data=data, method="REML")#specifies that residuals are correlated by site with compound symmetry
    summary(g4)
    g5<-lme(size.end~size.start+light, random=~1|site, na.action=na.omit, data=data, method="REML") #specifies site as a random effect
    summary(g5)
    anova(g4,g5)#models are identical
    anova(g3,g5)#split plot (light) not significant as random effect
    anova(g1,g5)#site as random effect not significant
    g6<-lme(size.end~size.start+light, random=list(~1|site, ~0+size.start|site), na.action=na.omit, data=data, method="REML")#intercept and effect of size varies by site, 
        #specifying in a list means R doesn't calculate correlation between random intercept and slope, which is its default (but not the default in SAS, I think)
    anova(g5,g6)#not a significant difference, size effect doesn't vary sig. by site
#Can include both random effecs and correlation in lme function, for example with repeated measures
    # e.g. correlation=corAR1(form=~year) would add auto-regressive order 1 correlation structure by year within random-effect groups ,if your data included multiple years
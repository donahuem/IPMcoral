rm(list=ls(all=TRUE)) ## Clear all
setwd("C://Users/Megan/Desktop/IPMMR184/latestripmstuff")# for MeganR
sfdat<-read.csv("EKBSFD.csv")
View(sfdat)
SFMC<-subset(sfdat,Species=="Mcap")
View(SFMC)
hist(SFMC$Ln.Area)

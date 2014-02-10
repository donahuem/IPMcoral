library(pscl)# for zero inflated poisson
MCr$lnspcover<-log(1+MCr$spcover_t)
rcrt.glm.SU <- glm(recruits~lnspcover*utrans,data=MCr,family=poisson)
drop1(rcrt.glm.SU,test="Chisq")
#drop spcover:utrans p-value=0.09382
rcrt.glm.SU0 <- glm(recruits~lnspcover+utrans,data=MCr,family=poisson)
drop1(rcrt.glm.SU0,test="Chisq")
#keep lnspcover p-value=1.203e-15 and utrans p-value <2.2e-16
summary(rcrt.glm.SU0)
#overdispersed Resdeviance/df=403.47/229=1.762


##Tried quasipoisson and it does not improve overdispersion (identical)
rcrt.glm.SU <- glm(recruits~lnspcover*utrans,data=MCr,family=quasipoisson)
summary(rcrt.glm.SU)


##Tried negative binomial
rcrt.nb.SU <- glm.nb(recruits~lnspcover*utrans,data=MCr)

drop1(rcrt.nb.SU,test="Chisq")
#drop lnspcover:utrans p-value=0.4408
rcrt.nb.SU0 <- glm.nb(recruits~lnspcover+utrans,data=MCr)
drop1(rcrt.nb.SU0,test="Chisq")
#keep lnspcover p-value=1.973e-09 and utrans p-value=4.431e-12
summary(rcrt.nb.SU0)
#does decrease residual deviance/df=195.52/229=0.854


###Try zero inflated models
rcrt.zip.SU0 <- zeroinfl(recruits~lnspcover+utrans,data=MCr,dist="poisson", link="logit")
rcrt.glm.SU0 <- glm(recruits~lnspcover+utrans,data=MCr,family=poisson)

AIC(rcrt.zip.SU0,rcrt.glm.SU0)
##zero inflated is a better fit by AIC delta=68.2156

rcrt.zinb.SU0 <- zeroinfl(recruits~lnspcover+utrans,data=MCr,dist="negbin", link="logit")
rcrt.nb.SU0 <- glm.nb(recruits~lnspcover+utrans,data=MCr)
AIC(rcrt.zinb.SU0,rcrt.nb.SU0)
##zero inflated is a better fit by AIC delta=1.3636

##negative binomial is a better fit in zero inflated and non inflated models



par(mfrow=c(2,2))
plot(resid(rcrt.glm.SU0)~MCr$utrans,ylim=c(-2,8))
plot(resid(rcrt.nb.SU0)~MCr$utrans,ylim=c(-2,8))
plot(resid(rcrt.zip.SU0)~MCr$utrans,ylim=c(-2,8))
plot(resid(rcrt.zinb.SU0)~MCr$utrans,ylim=c(-2,8))


xyplot(resid(rcrt.glm.SU0)~MCr$utrans,ylim=c(-2,8))
xyplot(resid(rcrt.nb.SU0)~MCr$utrans,ylim=c(-2,8))
xyplot(resid(rcrt.zip.SU0)~MCr$utrans,ylim=c(-2,8))
xyplot(resid(rcrt.zinb.SU0)~MCr$utrans,ylim=c(-2,8))

##Residuals look best for nb with lnspcover and utrans main effects
par(mfrow=c(2,2))
plot(MCr$lnspcover,MCr$recruits,col=cl[MCr$utrans],pch=18,main="poisson")
for (i in 1:10){
  lines(y,exp(coefficients(rcrt.glm.SU0)[1]+coefficients(rcrt.glm.SU0)[(1+i)]*y),col=cl[i])
}
plot(MCr$lnspcover,MCr$recruits,col=cl[MCr$utrans],pch=18,main="negative binomial")
for (i in 1:10){
  lines(y,exp(coefficients(rcrt.nb.SU0)[1]+coefficients(rcrt.nb.SU0)[(1+i)]*y),col=cl[i])
}
plot(MCr$lnspcover,MCr$recruits,col=cl[MCr$utrans],pch=18,main="zero-inflated poisson")
for (i in 1:10){
  lines(y,exp(coefficients(rcrt.zinb.SU0)[1]+coefficients(rcrt.zinb.SU0)[(1+i)]*y),col=cl[i])
}
plot(MCr$lnspcover,MCr$recruits,col=cl[MCr$utrans],pch=18,main="zero-inflated negative binomial")
for (i in 1:10){
lines(y,exp(coefficients(rcrt.zinb.SU0)[1]+coefficients(rcrt.zinb.SU0)[(1+i)]*y),col=cl[i])
}
legend("topleft", legend = sprintf(levels(MCr$utrans)), col=cl, lty=1,lwd=2, xjust = 1, bg = "white")
r1<-c(rep(0,10))
for (i in 1:10){
  r1[i]<-exp(coefficients(rcrt.nb.SU0)[1]+coefficients(rcrt.nb.SU0)[(1+i)]*y)
}
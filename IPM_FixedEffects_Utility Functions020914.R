#################################MeganR Fixed Effect Utility Functions#############################################################################

## utility functions
#probability of surviving
s.x=function(x,params) {
  u=exp(params[1,site]+params[2,site]*x)
  return(u/(1+u))
}

#growth functions by fateIPM
gg.yx=function(xp,x,params) {
  g<-dnorm(xp,params[3,site]+params[4,site]*x+params[5,site]*exp((-1)*params[6,site]*x),params[7,site])
}

#extreme shrinkage growth
gs.yx=function(xp,x,params) {
  g<-dnorm(xp,params[8,site]+params[9,site]*x,params[10,site])
}

ig.yx=function(xp,x,params) {
  g<-dnorm(xp,params[11,site]+params[12,site]*x,params[13,site])
}

ug.yx=function(xp,x,params) {
  g<-dnorm(xp,params[14,site]+params[15,site]*x,params[16,site])
}

# probability of extreme shrinkage    
p.es.x=function(x,params) {
  u=exp(params[17,site])
  return(u/(1+u))
}

# probability of fissing     
p.fiss.x=function(x,params) {
  u=exp(params[18,site]+params[19,site]*x)
  return(u/(1+u))
}
# probability of fusion
p.fuse.x=function(x,params) {
  u=exp(params[20,site]+params[21,site]*x)
  return(u/(1+u))
}
#probability of large fusion  
p.fuseLG.x=function(x,params) {
  u=exp(params[22,site]+params[23,site]*x)
  return(u/(1+u))
}

#Fission offspring using separate vglm pospoisson function for hi vs lo number of fission products
c.yx=function(xp,x,params) {
  nfiss<-exp(params[25,site]+params[26,site]*x)
  ##return: prob of fission * 1-prob of lo number*number of fission products * prob that a fisser of size x prodcues a product of size xp
  p.fiss.x(x,params)*nfiss*dunif(xp,min=-1.908971,max=x)
}

# number of clones per adult for "constant correction" of clonal matrix (no size distribution of offspring)
c.y=function(x,params) {
  nfiss<-exp(params[25,site]+params[26,site]*x)
  u=p.fiss.x(x,params)*nfiss
  return (u)
}

#recruitment = sexual offspring 

f.yx=function(xp,x,params){
  params[29,site]*x*dnorm(xp,params[27,site],params[28,site])
}

# recruits per adult for "constant correction" of fecundity matrix (no size distribution of offspring)
f.y=function(xp,params) {
  params[29,site]*xp
  
}

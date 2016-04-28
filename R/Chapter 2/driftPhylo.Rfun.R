########################################
#simulating phylogenetic data###########
########################################
require(phytools)
require(plyr); require(dplyr)
require(psych)
require(doParallel)
require(mvtnorm)
#registerDoParallel(cores = round(detectCores()/2))

# W- Pooled within covariance matrix
# means- species means
# phy- phylogeny
# n.s- vector of sample sizes
# inter- number of interaction
# parallel- run function on parallel (recomended)

sim_phylo<-function(vcv, means, phy, n.s, inter=1000, parallel=FALSE){
  pics<-apply(means, 2,function(x) pic(x, phy))
  B.pic<-var(pics)
  B<-var(means)
  n<-length(phy$tip.label)
  
  bm.sim<-function(z) {
    data.s<-sim.corrs(phy, vcv*tr(B.pic)/tr(vcv))
    ws<-rmvnorm(sum(n.s),sigma=vcv)
    sps<-rep(rownames(data.s),times=n.s)
    de<-data.frame(sps=factor(sps,unique(sps)),ws) %>% 
      group_by(.,sps) %>% summarize_each(funs(mean))
    data.s<-data.s+de[,-1]
    list(data.s, ws)
  }
  
  data.s<-alply(1:inter, 1, bm.sim, .parallel = parallel)
  list(sims=data.s,vcv=vcv,phy=phy,n=n, n.s=n.s, pics=pics,means=means)
}

#sim- output from the function sim_phylo 
power_phylo<-function(sim, parallel=FALSE){
  data.s<-sim$sims
  W<-sim$vcv
  n<-sim$n
  phy<-sim$phy
  eigenW<-eigen(W)
  inter<-length(data.s)
  
  typeI.WB<-ldply(data.s, function(x1){
    w<-var(x1[[2]])
    b<-var(x1[[1]])
    
    eigenw<-eigen(w)
    pcBw.s<-as.matrix(x1[[1]]) %*% eigenw$vectors
    v.pcBw.s<-var(pcBw.s)
    df.s  <- data.frame(W=log(eigenw$value),
                        B=log(diag(v.pcBw.s)))
    sltest<-lm(B~W, data=df.s)%>%
      linearHypothesis(., "W=1")
    sltest<-sltest[[6]][2]<c(0.05,0.01)
    
    C<-cov2cor(v.pcBw.s)
    dimnames(C)<-list(1:dim(C)[1],1:dim(C)[1])
    ctest<-corr.p(r= C, n = n, adjust = "bonferroni")
    ctest05<-sum(ctest$p[upper.tri(ctest$p)]<0.05)
    ctest01<-sum(ctest$p[upper.tri(ctest$p)]<0.01)
    
    sumcors<-C[lower.tri(C)] %>% abs %>% sort
    sumcors<-data.frame(mean=mean(sumcors), sd=sd(sumcors))
    
    if(dim(x1[[2]])[1]<dim(x1[[2]])[2]) C<-C[1:n,1:n]
    
    ctest.r<-corr.p(r= C, n = n, adjust = "bonferroni")
    ctest05.r<-sum(ctest.r$p[upper.tri(ctest.r$p)]<0.05)
    ctest01.r<-sum(ctest.r$p[upper.tri(ctest.r$p)]<0.01)
    
    data.frame(slt05=sltest[1], 
               slt01=sltest[2],
               cor05=as.numeric(ctest05),
               cor05.r=as.numeric(ctest05.r),
               cor01=ctest01,
               cor01.r=ctest01.r,
               cormean=sumcors[1],
               corsd  =sumcors[2])
  },.parallel = parallel)
  
  typeI.B<-ldply(data.s, function(x1){
    w<-W
    b<-var(x1[[1]])
    
    eigenw<-eigen(w)
    pcBw.s<-as.matrix(x1[[1]]) %*% eigenw$vectors
    v.pcBw.s<-var(pcBw.s)
    df.s  <- data.frame(W=log(eigenw$value),
                        B=log(diag(v.pcBw.s)))
    sltest<-lm(B~W, data=df.s)%>%
      linearHypothesis(., "W=1")
    sltest<-sltest[[6]][2]<c(0.05,0.01)
    
    C<-cov2cor(v.pcBw.s)
    dimnames(C)<-list(1:dim(C)[1],1:dim(C)[1])
    ctest<-corr.p(r= C, n = n, adjust = "bonferroni")
    ctest05<-sum(ctest$p[upper.tri(ctest$p)]<0.05)
    ctest01<-sum(ctest$p[upper.tri(ctest$p)]<0.01)
    
    sumcors<-C[lower.tri(C)] %>% abs %>% sort
    sumcors<-data.frame(mean=mean(sumcors), sd=sd(sumcors))
    
    if(dim(x1[[2]])[1]<dim(x1[[2]])[2]) C<-C[1:n,1:n]
    
    ctest.r<-corr.p(r= C, n = n, adjust = "bonferroni")
    ctest05.r<-sum(ctest.r$p[upper.tri(ctest.r$p)]<0.05)
    ctest01.r<-sum(ctest.r$p[upper.tri(ctest.r$p)]<0.01)
    
    data.frame(slt05=sltest[1], 
               slt01=sltest[2],
               cor05=as.numeric(ctest05),
               cor05.r=as.numeric(ctest05.r),
               cor01=ctest01,
               cor01.r=ctest01.r,
               cormean=sumcors[1],
               corsd  =sumcors[2])
  },.parallel = TRUE)
  
  typeI.WB.pic<-ldply(data.s, function(x1){
    w<-var(x1[[2]])
    pic<-apply(x1[[1]],2,function(a) pic(a,phy))
    b<-var(pic)
    
    eigenw<-eigen(w)
    
    pcBw.s<-pic %*% eigenw$vectors
    v.pcBw.s<-var(pcBw.s)
    df.s  <- data.frame(W=log(eigenW$value),
                        B=log(diag(v.pcBw.s)))
    sltest<-lm(B~W, data=df.s)%>%
      linearHypothesis(., "W=1")
    sltest05<-sltest[[6]][2]<0.05
    sltest01<-sltest[[6]][2]<0.01
    
    C<-cov2cor(v.pcBw.s)
    dimnames(C)<-list(1:dim(C)[1],1:dim(C)[1])
    ctest<-corr.p(r= C, n = n-1, adjust = "bonferroni")
    ctest05<-sum(ctest$p[upper.tri(ctest$p)]<0.05)
    ctest01<-sum(ctest$p[upper.tri(ctest$p)]<0.01)
    
    sumcors<-C[lower.tri(C)] %>% abs %>% sort
    sumcors<-data.frame(mean=mean(sumcors), sd=sd(sumcors))
    
    if(dim(pic)[1]<dim(x1[[1]])[2]) C<-C[1:(n-1),1:(n-1)]
    
    ctest.r<-corr.p(r= C, n = n-1, adjust = "bonferroni")
    ctest05.r<-sum(ctest.r$p[upper.tri(ctest.r$p)]<0.05)
    ctest01.r<-sum(ctest.r$p[upper.tri(ctest.r$p)]<0.01)
    
    data.frame(slt05=sltest05,
               slt01=sltest01,
               cor05=as.numeric(ctest05),
               cor05.r=as.numeric(ctest05.r),
               cor01=ctest01,
               cor01.r=ctest01.r,
               cormean=sumcors[1],
               corsd  =sumcors[2])
  },.parallel = parallel)
  
  typeI.B.pic<-ldply(data.s, function(x1){
    w<-W
    pic<-apply(x1[[1]],2,function(a) pic(a,phy))
    b<-var(pic)
    
    eigenw<-eigen(w)
    
    pcBw.s<-pic %*% eigenw$vectors
    v.pcBw.s<-var(pcBw.s)
    df.s  <- data.frame(W=log(eigenW$value),
                        B=log(diag(v.pcBw.s)))
    sltest<-lm(B~W, data=df.s)%>%
      linearHypothesis(., "W=1")
    sltest05<-sltest[[6]][2]<0.05
    sltest01<-sltest[[6]][2]<0.01
    
    C<-cov2cor(v.pcBw.s)
    dimnames(C)<-list(1:dim(C)[1],1:dim(C)[1])
    ctest<-corr.p(r= C, n = n-1, adjust = "bonferroni")
    ctest05<-sum(ctest$p[upper.tri(ctest$p)]<0.05)
    ctest01<-sum(ctest$p[upper.tri(ctest$p)]<0.01)
    
    sumcors<-C[lower.tri(C)] %>% abs %>% sort
    sumcors<-data.frame(mean=mean(sumcors), sd=sd(sumcors))
    
    if(dim(pic)[1]<dim(x1[[1]])[2]) C<-C[1:(n-1),1:(n-1)]
    
    ctest.r<-corr.p(r= C, n = n-1, adjust = "bonferroni")
    ctest05.r<-sum(ctest.r$p[upper.tri(ctest.r$p)]<0.05)
    ctest01.r<-sum(ctest.r$p[upper.tri(ctest.r$p)]<0.01)
    
    data.frame(slt05=sltest05,
               slt01=sltest01,
               cor05=as.numeric(ctest05),
               cor05.r=as.numeric(ctest05.r),
               cor01=ctest01,
               cor01.r=ctest01.r,
               cormean=sumcors[1],
               corsd  =sumcors[2])
  },.parallel = parallel)
  
  TYPEI.WB<-
    data.frame(error="WB",
               n,
               slError05=sum(typeI.WB$slt05)/inter,
               crError05=sum(typeI.WB$cor05>0)/inter,
               crError05.r=sum(typeI.WB$cor05.r>0)/inter,
               slError01=sum(typeI.WB$slt01)/inter,
               crError01=sum(typeI.WB$cor01>0)/inter,
               crError01.r=sum(typeI.WB$cor01.r>0)/inter,
               crMean=mean(typeI.WB$mean),
               crSD=mean(typeI.WB$sd),
               pic.slError05=sum(typeI.WB.pic$slt05)/inter,
               pic.crError05=sum(typeI.WB.pic$cor05>0)/inter,
               pic.crError05.r=sum(typeI.WB.pic$cor05.r>0)/inter,
               pic.slError01=sum(typeI.WB.pic$slt01)/inter,
               pic.crError01=sum(typeI.WB.pic$cor01>0)/inter,
               pic.crError01.r=sum(typeI.WB.pic$cor01.r>0)/inter,
               pic.crMean=mean(typeI.WB.pic$mean),
               pic.crSD=mean(typeI.WB.pic$sd))
  
  TYPEI.B<-
    data.frame(error="B",
               n,
               slError05=sum(typeI.B$slt05)/inter,
               crError05=sum(typeI.B$cor05>0)/inter,
               crError05.r=sum(typeI.B$cor05.r>0)/inter,
               slError01=sum(typeI.B$slt01)/inter,
               crError01=sum(typeI.B$cor01>0)/inter,
               crError01.r=sum(typeI.B$cor01.r>0)/inter,
               crMean=mean(typeI.B$mean),
               crSD=mean(typeI.B$sd),
               pic.slError05=sum(typeI.B.pic$slt05)/inter,
               pic.crError05=sum(typeI.B.pic$cor05>0)/inter,
               pic.crError05.r=sum(typeI.B.pic$cor05.r>0)/inter,
               pic.slError01=sum(typeI.B.pic$slt01)/inter,
               pic.crError01=sum(typeI.B.pic$cor01>0)/inter,
               pic.crError01.r=sum(typeI.B.pic$cor01.r>0)/inter,
               pic.crMean=mean(typeI.B.pic$mean),
               pic.crSD=mean(typeI.B.pic$sd))
  
  return(rbind(TYPEI.B,TYPEI.WB))
}

#data frame contains:
#error - which matrices were resampled. B- only B matrix, W considered fixed; WB- both W and B matrices were resampled.
#n - species sample size
#slError_ - Slope type I error
#crError_ - Correlation type I error
#_Error05 - error rate evaluated for alpha=0.05
#_Error01 - error rate evaluated for alpha=0.01
#crError_.r - error rate for reduced number of PCs (for n<number of variables)
#cr.Mean,crSD - mean and standard deviation of the absolute value of simulated correlations
#pic._Error_ - error estimated on Phylogenetic Independent Contrasts

simTests_phylo<-function(sim, parallel=FALSE, pic=FALSE){
  data.s<-sim$sims
  n<-sim$n
  phy<-sim$phy
  inter<-length(data.s)
  tests<-llply(data.s, function(x1){
    w<-var(x1[[2]])
    
    if(!pic) {
      b<-var(x1[[1]])
      } else {
        pic<-apply(x1[[1]],2,function(a) pic(a,phy))
        b<-var(pic)
        }
    
    eigenw<-eigen(w)
    pcBw.s<-as.matrix(x1[[1]]) %*% eigenw$vectors
    v.pcBw.s<-var(pcBw.s)
    df.s  <- data.frame(W=log(eigenw$value),
                        B=log(diag(v.pcBw.s)))
    sltest<-lm(B~W, data=df.s)

    C<-cov2cor(v.pcBw.s)

    return(list(sltest=sltest,C=C))
    },.parallel = parallel)
  sltests<-llply(tests, function(x) x$sltest)
  Cs<-llply(tests, function(x) x$C)
  return(list(sltests=sltests,Cs=Cs))
}

# a<-FAMILIES$Mustelidae
# sims<-sim_phylo(vcv = a$G,
#                 means = a$medias,
#                 phy = a$phy,
#                 n.s = laply(a$main,function(z) z$lmodel$df.residual)+1,
#                 inter = 200,
#                 parallel = TRUE)
# 
# simTests<-simTests_phylo(sims, parallel = FALSE, pic=TRUE)
# 
# #slope empirical confint
# simTests$sltests %>% 
#   laply(., function(x) coef(x)["W"]) %>% 
#   sort %>% 
#   quantile(.,probs = c(0.025,0.975))
# 
# #slope semi-empirical confint
# m<-simTests$sltests %>% 
#   laply(., function(x) coef(x)["W"]) %>% mean
# s<-simTests$sltests %>% 
#   laply(., function(x) coef(x)["W"]) %>% sd
# li<-qnorm(0.025) #inferior quantile
# ui<-qnorm(0.975) #upper quantile
# data.frame("2.5%"=m+li*s, "97.5%"=m+ui*s)
# 
# #correlation empirical threshold 
# Cs<-simTests$Cs %>% laply(., function(x) abs(x)) #put absolute correlations in an array
# Cs<-apply(Cs,2:3, sort)
# Cs[0.95*200,,]
# mean(Cs[0.95*200,,][upper.tri(diag(35))]) #
# 
# #correlation semi-empirical threshold 
# Csm <-apply(Cs,2:3, mean) #mean absolute correlations
# Cssd<-apply(Cs,2:3, sd)   #sd absolute correlations
# ui<-qnorm(0.95) #upper quantile
# Csm+Cssd*ui
# mean((Csm+Cssd*ui)[upper.tri(diag(35))])

#### from Proa et al, Evolution, 2013:
# TYPE I ERROR RATES FOR TESTING GENETIC DRIFT WITH PHENOTYPIC COVARIANCE MATRICES: A SIMULATION STUDY

simulateGP<-function(method, m, tNe, pop, n, sim.n, G, P)
  {
  #method = choose "sim1", "sim2", "sim3", "sim4", "sim5", according to description in manuscript and bellow
  #m = genetic covariance matrix dimension
  #tNe = number of generations (t) divided by effective population size (Ne)
  #pop = number of populations
  #n = sample size for each population
  #sim.n = number of simulations

    #define an empty variable to contain the results
  beta <- vector()
  prob <- vector()
  require(MASS)
  require(clusterGeneration)
  
   for (i in 1:sim.n) 
  {
         #if (method=="sim1")
         #Predetermined G and P will be used
         if (method=="sim2") {  P <- runif(1,min=1,max=10)*G   }
    #generate samples from pop populations with n obs each. the ancestral vector is composed of m zeros.
    M <-mvrnorm(pop, rep(0, ncol(G)), tNe*G)
    group <- factor( rep( seq(1:pop), each = n))
    data <- matrix(0, pop*n, ncol(P) )
    
    for (j in 1: (n*pop) ) {data[j,] <- mvrnorm(1, M[group[j],], P)}

    data<-adply(1: (n*pop), 1, function(j) mvrnorm(1, M[group[j],], P))[,-1]
    #calculate matrix of mean vectors from simulations
    
    zmeans <- matrix(unlist(by(data, group, colMeans)),  nrow=pop, ncol = ncol(G), byrow=T)
    
    #calculate within-group phenotypic covariance and extract eigenvalues
    
    eigW <- eigen(cov( data - zmeans[rep(1:nrow(zmeans), each = n),]))
    
    #Project mean vectors for the pop pops on within group eigenvectors
    
    zm.proj <- zmeans %*% eigW$vectors
    
    #calculate among-group variance
    
    v <- diag(cov(zm.proj))
    
    #perform Ackermann and Cheverud test with t-test for unity slope
    
    model <- summary(lm(log(v) ~ log(eigW$values)))
    prob[i] <- pt(abs( (model$coefficients[2,1]-1) / model$coefficients[2,2]), df = model$df[2], lower.tail=F)*2 # distribuição de student
    beta[i] <- model$coefficients[2,1]
  }
  #Calculate 95% confidence limits
  Beta <- data.frame(Lower = NA, Upper = NA, Mean = NA)
  Beta$Upper<-mean(beta)+1.96*sd(beta)
  Beta$Lower<-mean(beta)-1.96*sd(beta)
  Beta$Mean <- mean(beta)
    
  tIe <- mean(prob < 0.05)
  
  return(list("TypeIerror" = tIe,
              "Probabilities" = prob,
              "Beta.ic" = Beta,
              "betas.dist" = beta) )
}



simulateGPBeta<-function(method, m, tNe, pop, n, sim.n, G="NULL", P="NULL",smax){
  #define an empty variable to contain the results
  beta <- vector()
  require(MASS)
  require(clusterGeneration)
  
  for (i in 1:sim.n) 
  {
    if (method=="sim1")
    {
      if (G=="NULL") { stop("you need to provide a valid G matrix!\n") 
      }
      if (P=="NULL") {stop("you need to provide a valid P matrix!\n") 
      }
      #Predetermined G and P will be used
    }
    if (method=="sim2")
    {
      #G is defined as a random positive definite matrix, and P = xG, only different by a random constant of proportionality.
      #G<-genPositiveDefMat(dim=m,lambdaLow=1,ratioLambda=10)$Sigma
      P<-runif(1,min=1,max=10)*G
    }
    if (method=="sim3")
    {
      #G and P matrices are defined randomly and independently from each other, with the constraint that the variances in P are always larger than the respective variances in G.
      G<-genPositiveDefMat(dim=m,covMethod="unifcorrmat",rangeVar=c(1,10))$Sigma
      R<-genPositiveDefMat(dim=m,covMethod="unifcorrmat",rangeVar=c(1,1))$Sigma
      var<-runif(m,2,9)*diag(G)
      P<-diag(sqrt(var))%*%R%*%diag(sqrt(var))
    }
    if (method=="sim4")
    {
      #G and E matrices are defined randomly and independently from each other, whereas P = G+E. P and G do not share a common PC structure but are related
      G<-genPositiveDefMat(dim=m,covMethod="unifcorrmat",rangeVar=c(1,10))$Sigma
      P<-G+genPositiveDefMat(dim=m,covMethod="unifcorrmat",rangeVar=c(0,smax))$Sigma
    }
    if (method=="sim5")
    {
      #G and E matrices are defined randomly and independently from each other, whereas P = G+E. P and G share a common PC structure
      G<-genPositiveDefMat(dim=m,lambdaLow=1,ratioLambda=10)$Sigma
      P<-G+genPositiveDefMat(dim=m,covMethod="unifcorrmat",rangeVar=c(0,smax))$Sigma
    }
    
    #generate samples from pop populations with n obs each. the ancestral vector is composed of m zeros.
    M <- mvrnorm(pop, rep(0, ncol(G)), tNe*G)
    group<-factor(rep( seq(1:pop), each=n))
    data <- matrix(0, pop*n, ncol(P))
    for (j in 1:(n*pop)) {data[j,] <- mvrnorm(1, M[group[j],], P)}
    
    #calculate matrix of mean vectors from simulations
    
    zmeans <- matrix(unlist(by(data,group,colMeans)),nrow=pop,ncol=ncol(G),byrow=T)
    
    #calculate within-group phenotypic covariance and extract eigenvalues
    
    eigW <- eigen(cov(data-zmeans[rep(1:nrow(zmeans),each=n),]))
    
    #Project mean vectors for the 15 pops on within group eigenvectors
    
    zm.proj <- zmeans %*% eigW$vectors
    
    #calculate among-group variance
    
    v <- diag(cov(zm.proj))
    
    #calculate slope of Ackermann and Cheverud test
    
    model<-summary(lm(log(v)~log(eigW$values)))
    beta[i]<-model$coefficients[2,1]
  }
  #Calculate 95% confidence limits
  UL<-mean(beta)+1.96*sd(beta)
  LL<-mean(beta)-1.96*sd(beta)
  return(c(mean(beta),UL,LL))
}

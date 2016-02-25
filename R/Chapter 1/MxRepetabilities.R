#Montando as tabelas de Repetibilidades


# faltou calcular as repetibilidades por monte carlo
AddMCReps <- function(x){
  x[[16]][[2]] <- matrix(NA, nrow = 1, ncol= 2)
  x[[16]][[2]][1] <- if(x[[9]]>14) MonteCarloRep( x[[11]][[1]], sample.size = x[[9]], ComparisonFunc = RandomSkewers, iterations = 1000, parallel = T) else NA
  x[[16]][[2]][2] <- if(x[[9]]>14) MonteCarloRep( x[[11]][[1]], sample.size = x[[9]], ComparisonFunc = KrzCor, iterations = 1000, parallel = T) else NA
  #x[[16]][[2]][3] <- if(x[[9]]>14) MonteCarloRep( x[[11]][[1]], sample.size = x[[9]], ComparisonFunc = PCAsimilarity, iterations = 1000, parallel = T) 
  #x[[16]][[2]][4] <- if(x[[9]]>14) MonteCarloRep( x[[11]][[1]], sample.size = x[[9]], ComparisonFunc = MantelCor, iterations = 1000, correlation = TRUE, parallel = T) 
  #x[[16]][[2]][5] <- if(x[[9]]>14) MonteCarloRep( x[[11]][[1]], sample.size = x[[9]], ComparisonFunc = KrzCor, iterations = 1000, correlation = TRUE, parallel = T)
  colnames(x[[16]][[2]])[1:2] <- c('rs', 'krz')
  names(x[[16]])[2] <- c('MCRep')
  return(x)
  
}
registerDoParallel(cores = 2)
current.data <- sp.main.data
sp.main.data.t <- current.data %>% llply(AddMCReps, .progress = "text")

mx.rep <- current.data %>% ldply(function(x) x$Mx.Rep$BootsRep) 
names(mx.rep) <- c("Especie", "BS.rs", "BS.krz", "pcas", "cor.mantel", "cor.krz")
mx.rep.mc <- sp.main.data.t%>% ldply(function(x) x$Mx.Rep$MCRep) 
names(mx.rep.mc) <- c("Especie", "MC.rs", "MC.krz")
tab.rep.mx <- cbind(mx.rep[,1:3], mx.rep.mc[,2:3])
tab.rep.mx$Especie %<>% gsub("_", ' ', .)
row.names(tab.rep.mx) <- tab.rep.mx[,1]
tab.rep.mx$n <- n.size$V1

cov.mx <- current.data %>% llply(function(x) x$matrix$cov)
mask <- ldply(cov.mx, function(x) !is.na(x[1]))[,2]

rep.list <- mx.rep[mask, ]
row.names(rep.list) <- rep.list$Especie
rep.list$Especie <- as.character(rep.list$Especie)
rep.list[dim(rep.list)[1] +1, ] <- c("Saguinus_P", rep(0.97, 5) )
rep.list[dim(rep.list)[1] +1, ] <- c("Saguinus_G", rep(0.75, 5) )
row.names(rep.list) <- rep.list$Especie
rep.list$Especie <- factor(rep.list$Especie, levels = unique(rep.list$Especie) )
rep.list$rs <- as.numeric(rep.list$rs)
rep.list$krz <- as.numeric(rep.list$krz)
rep.list$pcas <- as.numeric(rep.list$pcas)
rep.list$cor.mantel <- as.numeric(rep.list$cor.mantel)
rep.list$cor.krz <- as.numeric(rep.list$cor.krz)
str(rep.list)

####################### calculando o valor médio de comparação pra cada especie ##########

mean.sim<-function(x) {
  x[upper.tri(x)]<-t(x)[upper.tri(x)]
  diag(x)<-NA
  x<-rowMeans(x,na.rm=TRUE)
  return(x)
}
mean.compare.cov<- mx.compare[c(1:2)] %>% llply (function (x) x$correlations ) %>% ldply (function (x) mean.sim(x) )

dimnames(mean.compare.cov)[[1]] <- c(paste("Mean", mean.compare.cov[,1] ))
dimnames(mean.compare.cov)[[2]] <- dimnames(mean.compare.cov)[[2]]

row.names(mean.compare.cov)

mean.compare.cov<- t(mean.compare.cov[,-1])
tab.rep.mx.mix<- cbind(tab.rep.mx[mask,-1] , mean.compare.cov[-c(43:44),] )

xtable(tab.rep.mx.mix, digits = 3, caption = "Matrix repetabilities calculated via Bootstrap and Monte Carlo using RS and KRZ methods of comparisson.")




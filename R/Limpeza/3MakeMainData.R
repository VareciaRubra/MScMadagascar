################ ?Organizando tudo em listas por ranking taxonomico ##############

current.data <- raw.main.data$Eulemur_albifrons
makeMainData <- function (current.data, final = FALSE) {
  x = vector("list", 13)
  x[[1]] <- dplyr::select(current.data, c(Arquivo:Take)) # pegando as info de cada replica
  x[[2]] <- dplyr::select(current.data, c(IS_PM:BA_OPI)) # pegando as ed de cada replica
  x[[3]] <- base::unique(select(current.data, c(Arquivo:Data_dado))) # pegando as info únicas de cada indivíduo
  x[[4]] <- ddply(dplyr::select(current.data, c(Tombo, IS_PM:BA_OPI)), .(Tombo), numcolwise(mean)) # calculando as médias individuais, excluindo a primeira coluna que ta a info de id
  rownames(x[[4]]) <- x[[4]][,1] # nomeando as linhas como o numero de tombo para manter correspondencia
  x[[4]]<- x[[4]][, -1]
  names(x)[1:4] <- c('info.raw', 'ed.raw', 'info', 'ed') # nomemando as listas
  x[[5]] <- t(apply(x[[4]], 1, FUN = function(x) x/ exp( mean( as.matrix(log (x)) ) )  ) ) 
  rownames(x[[5]]) <- x[[3]]$Tombo # nomeando as linhas como o numero de tombo para manter correspondencia
  x[[6]] <- log(x[[4]])
  rownames(x[[6]]) <- rownames(x[[4]])
  x[[7]] <- t(apply(x[[4]], 1, FUN = function(x) exp( mean( as.matrix(log (x)) ) )  ) )  
  x[[8]] <- mean(x[[7]])
  names(x)[5:8] <- c('sizeless', 'log', 'gm.ind', 'gm.mean')
  x[[9]] <- dim(x[[3]])[1]
  x[[10]] <- colMeans(x[[4]])
  names(x)[9:10] <- c('sample.size', 'ed.means')
  x[[11]] <- vector("list", 4)
  x[[11]][[1]] <- if(x[[9]]>20) var(x[[4]]) else matrix(data = NA, nrow = 39, ncol = 39) 
  x[[11]][[2]] <- if(x[[9]]>20) cov2cor(x[[11]][[1]]) else matrix(data = NA, nrow = 39, ncol = 39) 
  x[[11]][[3]] <- if(x[[9]]>20) var(x[[5]]) else matrix(data = NA, nrow = 39, ncol = 39) 
  x[[11]][[4]] <- if(x[[9]]>20) var(x[[6]]) else matrix(data = NA, nrow = 39, ncol = 39) 
  names(x[[11]]) <- 'matrix'
  names(x[[11]])[1:4] <- c('cov','cor', 'cov.sizeless', 'cov.log')
if (compare.size == TRUE)
  {x[[12]] <- ldply(.data = x[[11]][-2], .fun = MeanMatrixStatistics, parallel = TRUE)
    names (x[[12]]) <- 'mean.mx.stats'
   x[[13]] <- vector("numeric", 4)
   x[[13]][[1]]<- if(x[[9]]>20) BootstrapRep( x[[4]], ComparisonFunc = RandomSkewers, iterations = 1000, parallel = T) else NA
   x[[13]][[2]]<- if(x[[9]]>20) BootstrapRep( x[[4]], ComparisonFunc = MantelCor, iterations = 1000, correlation = TRUE, parallel = T) else NA
   x[[13]][[3]]<- if(x[[9]]>20) BootstrapRep( x[[5]], ComparisonFunc = RandomSkewers, iterations = 1000, parallel = T) else NA
   x[[13]][[4]]<- if(x[[9]]>20) BootstrapRep( x[[6]], ComparisonFunc = RandomSkewers, iterations = 1000, parallel = T) else NA
   names(x)[13] <- c('BootsRep')
   #x[[13]] <- t(x[[13]])
   #dimnames(x[[13]]) <- c('cov.RS','cor.Mt', 'cov.sizeless.RS', 'cov.log.RS')
   x[[14]] <- RandomSkewers(cov.x = x[[11]][-2], num.vectors = 1000, repeat.vector = x[[13]][-2]) 
   names(x)[[14]] <- 'rs.size.comparisson'
   return(x)}
  if (final == TRUE)
  {x[[15]] <- vector("list", 5)
   x[[15]][[1]] <- if(x[[9]]>20) Rarefaction(x[[4]], RandomSkewers, iterations = 1000, num.reps = x[[9]], parallel = T) else NA
   x[[15]][[2]] <- if(x[[9]]>20) Rarefaction(x[[4]], KrzCor, iterations = 1000, num.reps = x[[9]], parallel = T) else NA
   x[[15]][[3]] <- if(x[[9]]>20) Rarefaction(x[[4]], PCAsimilarity, iterations = 1000, num.reps = x[[9]], parallel = T) else NA
   x[[15]][[4]] <- if(x[[9]]>20) Rarefaction(x[[4]], MatrixCor, correlation = TRUE, iterations = 1000, num.reps = x[[9]], parallel = T) else NA
   x[[15]][[5]] <- if(x[[9]]>20) Rarefaction(x[[4]], KrzCor, correlation = TRUE, iterations = 1000, num.reps = x[[9]], parallel = T) else NA
   names(x)[15] <- c('rarefaction')
   names(x[[15]])[1:5] <- c('rs', 'krz','pcas', 'cor.mantel', 'cor.krz')
   x[[16]] <- vector("list", 5)
   x[[16]][[1]] <- if(x[[9]]>20) BootstrapRep( x[[4]], ComparisonFunc = RandomSkewers, iterations = 1000, parallel = T) else NA
   x[[16]][[2]] <- if(x[[9]]>20) BootstrapRep( x[[4]], ComparisonFunc = KrzCor, iterations = 1000, parallel = T) else NA
   x[[16]][[3]] <- if(x[[9]]>20) BootstrapRep( x[[4]], ComparisonFunc = PCAsimilarity, iterations = 1000, parallel = T) else NA
   x[[16]][[4]] <- if(x[[9]]>20) BootstrapRep( x[[4]], ComparisonFunc = MantelCor, iterations = 1000, correlation = TRUE, parallel = T) else NA
   x[[16]][[5]] <- if(x[[9]]>20) BootstrapRep( x[[4]], ComparisonFunc = KrzCor, iterations = 1000, correlation = TRUE, parallel = T) else NA
   names(x)[16] <- c('BootsRep')
   names(x[[16]])[1:5] <- c('rs', 'krz','pcas', 'cor.mantel', 'cor.krz')
   return(x)}
  
  else return(x)
  
}

makeMainData(raw.main.data$Eulemur_albifrons, final = TRUE)

main.data <- llply(raw.main.data, makeMainData, final = TRUE, .progress = 'text')
all.main.data <- llply(all.raw.main.data,  final = FALSE, makeMainData)

Eulemur.image <- list (main.data = main.data,   
                         all.main.data = all.main.data)

save(Eulemur.image,
    file = arquivo.saida)



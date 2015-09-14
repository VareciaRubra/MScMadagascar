################ ?Organizando tudo em listas por ranking taxonomico ##############

current.data <- Sp.raw.main.data$Eulemur_rubriventer
makeMainData <- function (current.data, specie = TRUE, compare.size = FALSE, final = FALSE ) 
  {
  x = vector("list", 17 )
  x[1:11] <- NA
  x[[1]] <- dplyr::select(current.data, c(Arquivo:Take)) # pegando as info de cada replica
  x[[2]] <- dplyr::select(current.data, c(IS_PM:BA_OPI)) # pegando as ed de cada replica
  x[[3]] <- base::unique(select(current.data, c(Arquivo:Data_dado))) # pegando as info únicas de cada indivíduo
  x[[4]] <- ddply(dplyr::select(current.data, c(Tombo, IS_PM:BA_OPI)), .(Tombo), numcolwise(mean)) # calculando as médias individuais, excluindo a primeira coluna que ta a info de id
  rownames(x[[4]]) <- x[[4]][,1] # nomeando as linhas como o numero de tombo para manter correspondencia
  x[[4]]<- x[[4]][, -1]
  names(x)[1:4] <- c('info.raw', 'ed.raw', 'info', 'ed') # nomemando as listas
  x[[5]] <- t(apply(x[[4]], 1, FUN = function(x) x/ exp( mean( as.matrix(log (x)) ) )  ) ) 
  #rownames(x[[5]]) <- x[[3]]$Tombo # nomeando as linhas como o numero de tombo para manter correspondencia
  x[[6]] <- log(x[[4]])
  #rownames(x[[6]]) <- rownames(x[[4]])
  x[[7]] <- t(apply(x[[4]], 1, FUN = function(x) exp( mean( as.matrix(log (x)) ) )  ) )  
  x[[8]] <- mean(x[[7]])
  names(x)[5:8] <- c('sizeless', 'log', 'gm.ind', 'gm.mean')
  x[[9]] <- dim(x[[3]])[1]
  x[[10]] <- colMeans(x[[4]])
  names(x)[9:10] <- c('sample.size', 'ed.means')
  x[[11]] <- vector("list", 5)
    sp.number <- length(unique(x[[3]]$Especie) ) 
    sp.fit <-  if (specie == FALSE & sp.number > 1 & sum((table(x[[3]]$Especie)>=3) == FALSE) ==0  ) manova(as.matrix(x[[4]]) ~ Especie, data = as.data.frame(x[[3]])) else NA
    gr <- if(specie == TRUE | specie == FALSE & sp.number <= 1 | specie == FALSE & sum((table(x[[3]]$Especie)>=3) == FALSE) !=0    ) var(na.omit(x[[4]])) else CalculateMatrix(sp.fit) 
    pq <- matrix(data = NA, nrow = 39, ncol = 39) 
  x[[11]][[1]] <- if(x[[9]]>14) gr else pq
  x[[11]][[2]] <- if(x[[9]]>14) cov2cor(x[[11]][[1]]) else pq
  x[[11]][[3]] <- if(x[[9]]>14) var(na.omit(x[[5]]) ) else pq
  x[[11]][[4]] <- if(x[[9]]>14) var(na.omit(x[[6]]) ) else pq
  x[[11]][[5]] <- sp.fit 
  x[[11]][[6]] <- if( !is.na(x[[11]][[1]][1]) ) ExtendMatrix(cov.matrix = x[[11]][[1]], var.cut.off = 1e-04, ret.dim = NULL) else pq
  names(x)[11] <- 'matrix'
  names(x[[11]])[1:6] <- c('cov','cor', 'cov.sizeless', 'cov.log', 'fit', 'ext.cov.mx')

  x[[13]] <- vector("list", 4)
  names(x[[13]])[1:4] <-  c('cov','cor', 'cov.sizeless', 'cov.log')
  if (compare.size == TRUE)
  {x[[12]] <- if(x[[9]]>40) ldply(.data = x[[11]][-c(2,5,6)], .fun = MeanMatrixStatistics, parallel = TRUE, full.results = FALSE) else NA
   x[[13]][1]<- if(x[[9]]>14) BootstrapRep( x[[4]], ComparisonFunc = RandomSkewers, iterations = 1000, parallel = TRUE) else NA
   x[[13]][2]<- if(x[[9]]>14) BootstrapRep( x[[4]], ComparisonFunc = MantelCor, iterations = 1000, correlation = TRUE, parallel = T) else NA
   x[[13]][3]<- if(x[[9]]>14) BootstrapRep( x[[5]], ComparisonFunc = RandomSkewers, iterations = 1000, parallel = TRUE) else NA
   x[[13]][4]<- if(x[[9]]>14) BootstrapRep( x[[6]], ComparisonFunc = RandomSkewers, iterations = 1000, parallel = TRUE) else NA
   x[[14]] <- if(x[[9]]>14) RandomSkewers(cov.x = x[[11]][c(1,3,4)], num.vectors = 1000, repeat.vector = unlist(x[[13]])[c(1,3,4)]  ) else NA
  }   
  names (x)[12:14] <- c('mean.mx.stats', 'size.mx.rep.rs', 'rs.size.comparisson')
  
  x[[15]] <- vector("list", 5)
  x[[16]] <- vector("list", 1)
  x[[16]][[1]] <- matrix(NA, nrow = 1, ncol= 5)
  x[[17]] <- vector ("list", 3)
  if (final == TRUE)
  {
    x[[15]][[1]] <- if(x[[9]]>14) Rarefaction(x[[4]], RandomSkewers, iterations = 1000, num.reps = x[[9]], parallel = T) else NA
    x[[15]][[2]] <- if(x[[9]]>14) Rarefaction(x[[4]], KrzCor, iterations = 1000, num.reps = x[[9]], parallel = T) else NA
    x[[15]][[3]] <- if(x[[9]]>14) Rarefaction(x[[4]], PCAsimilarity, iterations = 1000, num.reps = x[[9]], parallel = T) else NA
    x[[15]][[4]] <- if(x[[9]]>14) Rarefaction(x[[4]], MatrixCor, correlation = TRUE, iterations = 1000, num.reps = x[[9]], parallel = T) else NA
    x[[15]][[5]] <- if(x[[9]]>14) Rarefaction(x[[4]], KrzCor, correlation = TRUE, iterations = 1000, num.reps = x[[9]], parallel = T) else NA
        
    x[[16]][[1]][1] <- if(x[[9]]>14) BootstrapRep( x[[4]], ComparisonFunc = RandomSkewers, iterations = 1000, parallel = T) else NA
    x[[16]][[1]][2] <- if(x[[9]]>14) BootstrapRep( x[[4]], ComparisonFunc = KrzCor, iterations = 1000, parallel = T) else NA
    x[[16]][[1]][3] <- if(x[[9]]>14) BootstrapRep( x[[4]], ComparisonFunc = PCAsimilarity, iterations = 1000, parallel = T) else NA
    x[[16]][[1]][4] <- if(x[[9]]>14) BootstrapRep( x[[4]], ComparisonFunc = MantelCor, iterations = 1000, correlation = TRUE, parallel = T) else NA
    x[[16]][[1]][5] <- if(x[[9]]>14) BootstrapRep( x[[4]], ComparisonFunc = KrzCor, iterations = 1000, correlation = TRUE, parallel = T) else NA
    
    x[[17]][[1]] <- if(x[[9]]>14) BootstrapR2(ind.data= x[[4]], iterations = 1000, parallel = TRUE) else NA
    x[[17]][[2]] <- if(x[[9]]>14) BootstrapR2(ind.data= x[[5]], iterations = 1000, parallel = TRUE) else NA
    x[[17]][[3]] <- if(x[[9]]>14) BootstrapR2(ind.data= x[[6]], iterations = 1000, parallel = TRUE) else NA
    
  }  
  
  names(x)[15] <- c('rarefaction')
  names(x[[15]])[1:5] <- c('rs', 'krz','pcas', 'cor.mantel', 'cor.krz')
  names(x)[16] <- c('Mx.Rep')
  names(x[[16]])[[1]] <- c('BootsRep')
  colnames(x[[16]][[1]])[1:5] <- c('rs', 'krz','pcas', 'cor.mantel', 'cor.krz')
  names(x)[17] <- 'BootsR2'
  names(x[[17]])[1:3] <- c('ed.means', 'sizeless','log')
      
  return(x)
  
}





#############################################
############# PRIMEIRA PARTE ################
####### OLHANDO PRO PLOT DE PC1 X PC2 #######
#############################################



extant.main.data<- raw.data %>% filter( Status != "Extinct") %>% filter( Genero != "Tarsius")
extant.main.data<- dlply(extant.main.data, .(All), tbl_df)
extant.main.data<- llply(extant.main.data, specie = TRUE, final = FALSE, makeMainData, .progress = 'text')

madagascar.main.data<- raw.data %>% filter( Regiao == "Madagascar") 
madagascar.main.data<- dlply(madagascar.main.data, .(All), tbl_df)
madagascar.main.data<- llply(madagascar.main.data, specie = TRUE, final = FALSE, makeMainData, .progress = 'text')

extant.madagascar.main.data<- raw.data %>% filter( Regiao == "Madagascar") %>% filter( Status != "Extinct") 
extant.madagascar.main.data<- dlply(extant.madagascar.main.data, .(All), tbl_df)
extant.madagascar.main.data<- llply(extant.madagascar.main.data, specie = TRUE, final = FALSE, makeMainData, .progress = 'text')

#all.main.data <- llply(All.raw.main.data, specie = TRUE, final = TRUE, makeMainData, .progress = 'text')
#sp.master.main.data <- llply(Sp.raw.main.data, specie = TRUE, final = FALSE, makeMainData, .progress = 'text', .inform = T)
#gen.master.main.data <- llply(Gen.raw.main.data, specie = FALSE, final = FALSE, makeMainData, .progress = 'text', .inform = T)


#############################################
############## SEGUNDA PARTE ################
#### ARQUIVO FINAL PARA SALVAR NO RData #####
#############################################

#makeMainData(raw.main.data[[1]], specie = TRUE, compare.size = TRUE, final = TRUE)

#main.data<- llply(raw.main.data, makeMainData, specie = TRUE, compare.size = FALSE, final = TRUE, .progress = 'text')

all.main.data<- llply(all.raw.main.data, makeMainData, specie = FALSE, compare.size = FALSE, final = TRUE, .progress = 'text')


#############################################
############# TERCEIRA PARTE ################
## Colocando Mx bem estimadas para os n<27###
#############################################
registerDoParallel(cores = 30)
sp.main.data <- llply(Sp.raw.main.data, specie = TRUE, compare.size = TRUE, final = TRUE, makeMainData, .progress = progress_text(char = "."), .inform = T)
gen.main.data <- llply(Gen.raw.main.data, specie = FALSE, compare.size = TRUE, final = TRUE, makeMainData, .progress = progress_text(char = "."), .inform = T)
  
Strepsirrhini.image.complete <- list (specie.lists = sp.main.data,   
                               genus.list = gen.main.data)
                             
save( Strepsirrhini.image.complete,
     file = "~/ataches/Strepsirrhini_image_complete_III.RData")

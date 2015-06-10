current.data <- raw.main.data$Hapalemur_griseus
DoubleRep <- function (current.data) {
  x = vector("list", 10)
  x[[1]] <- dplyr::select(current.data, c(Arquivo:Take)) # pegando as info de cada replica
  x[[2]] <- dplyr::select(current.data, c(IS_PM:BA_OPI)) # pegando as ed de cada replica
  x[[3]] <- unique(select(current.data, c(Arquivo:Data_dado))) # pegando as info únicas de cada indivíduo
  x[[4]] <- ddply(dplyr::select(current.data, c(Tombo, IS_PM:BA_OPI)), .(Tombo), numcolwise(mean))[,-1] # calculando as médias individuais, excluindo a primeira coluna que ta a info de id
  rownames(x[[4]]) <- x[[3]]$Tombo # nomeando as linhas como o numero de tombo para manter correspondencia
  x[[5]] <- ddply(dplyr::select(current.data, c(Tombo, IS_PM:BA_OPI)), .(Tombo), numcolwise(diff))[,-1]
  rownames(x[[5]]) <- x[[3]]$Tombo 
  names(x)[1:5] <- c('info.raw', 'ed.raw', 'info', 'ed', 'diff.ind') # nomemando as listas
  ###### tabelando quantas replicas diferentes tem por distancia (selecionando só os que sao T #################
  X2dists<- as.matrix(apply(x[[5]], 2, function (x)table(x!=0)))[2,]
  
  total <- dim(x[[3]])[1]
  quantos<- X2dists/total
  ####### matriz para indexar quais replicas sao diferentes entre si #############
  indexX2<- aggregate(x[[2]], by =list(x[[1]]$Tombo) , FUN = function (x) x[1] != x[2])
  dimnames(indexX2)[[1]]<- as.character(indexX2[,1])
  indexX2<- as.matrix(indexX2[,-1])
  X2distsdif<-as.matrix(apply(indexX2, 2, function (x)table(x == T)))[2,]
  cbind(X2dists, X2distsdif)
  
  for (i in 1:39 ) { 
    x2<-as.list(rep(colnames(x[[2]], X2dists[i][1])))
    names(x2)<- colnames(x[[2]])
    rm(i)
  }
  for (i in 1:39){
    x2[i]<- data.frame(rep(NA, X2dists[i][1])) 
    x2[[i]] <- as.factor(as.vector(dimnames(indexX2)[[1]][indexX2[, i]]))
    rm(i)
  }
  
  options(contrasts = c('contr.sum', 'contr.poly'), na.action = 'na.omit')
  clean.matrix <- x[[2]]
  for (i in 1:39)
  {
    clean.matrix [!(x[[1]]$Tombo %in% x2 [[i]]), i] <- NA
  }
  doub.rep<- CalcRepeatability(x[[1]]$Tombo, clean.matrix)
  all.rep<- CalcRepeatability(x[[1]]$Tombo, x[[2]])
  reps<- cbind(names(doub.rep) , as.data.frame(doub.rep), as.data.frame(all.rep), total,  X2dists, quantos  )
  
  return(reps)
}

mask.rep = laply(raw.main.data, count, Especie)[,2] > 4
rep.sp <- llply(raw.main.data[mask.rep], DoubleRep)
rep.gen <- ldply(all.raw.main.data, DoubleRep)
names(raw.main.data[[4]])


rep.sp %>% 
  ldply(function(x) x$doub.rep) %>% 
  {. %*% t(.)}

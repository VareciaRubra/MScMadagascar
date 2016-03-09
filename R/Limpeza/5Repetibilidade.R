current.data <- Sp.raw.main.data[[6]]
DoubleRep <- function (current.data, detailed = FALSE) {
  x = vector("list", 10)
  x[[1]] <- dplyr::select(current.data, c(Arquivo:Take)) # pegando as info de cada replica
  x[[2]] <- dplyr::select(current.data, c(IS_PM:BA_OPI)) # pegando as ed de cada replica
  x[[3]] <- base::unique(select(current.data, c(Arquivo:Data_dado))) # pegando as info únicas de cada indivíduo
  x[[4]] <- ddply(dplyr::select(current.data, c(Tombo, IS_PM:BA_OPI)), .(Tombo), numcolwise(mean)) # calculando as médias individuais, excluindo a primeira coluna que ta a info de id
  rownames(x[[4]]) <- x[[4]][,1] # nomeando as linhas como o numero de tombo para manter correspondencia
  x[[4]]<- x[[4]][, -1]
  x[[5]] <- ddply(dplyr::select(current.data, c(Tombo, IS_PM:BA_OPI)), .(Tombo), numcolwise( diff) ) 
  rownames(x[[5]]) <- x[[5]][,1] # nomeando as linhas como o numero de tombo para manter correspondencia
  x[[5]]<- x[[5]][, -1]
  names(x)[1:5] <- c('info.raw', 'ed.raw', 'info', 'ed', 'diff.ind') # nomemando as listas
  x[[6]] <- dim(x[[3]])[1]
  x[[7]] <- colMeans(x[[4]])
  
  names(x)[6:7] <- c('sample.size', 'ed.means')
  
  ###### tabelando quantas replicas diferentes tem por distancia (selecionando só os que sao T #################
  X2dists<- as.matrix(apply(x[[5]], 2, function (x) sum(x != 0)))
  total <- dim(x[[3]])[1]
  doub.rep<- X2dists
    
  ##### testes logicos pra ver se rola de calcular a repetibilidade:
  ###### 1. ver se tem diferença entre as replicas de pelo menos algum indivíduo, 
  ###### 2. ver se tem diferença entre as replicas de pelo menos algum indivíduo, 
  ###### 3. ver se tem mais de 4 bichos medidos com replica
  if( dim( table(x[[5]] !=0) ) !=1 & sum((x[[5]] !=0) == TRUE) !=0 & sum(X2dists) >= (39*4)  )
  { 
   ####### matriz para indexar quais replicas sao diferentes entre si #############
   indexX2<- aggregate(x[[2]], by =list(x[[1]]$Tombo) , FUN = function (x) x[1] != x[2])
   dimnames(indexX2)[[1]]<- as.character(indexX2[,1])
   indexX2<- as.matrix(indexX2[,-1])
   X2distsdif <- as.matrix(apply(indexX2, 2, function (x)sum(x == TRUE)))
   cbind(X2dists, X2distsdif)
   quantos<- X2distsdif/total
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
   if (detailed == TRUE) return(reps) else return(as.data.frame(doub.rep))
  } else return( as.data.frame(doub.rep)/NA )
   
} 

mask.rep = ldply(Sp.raw.main.data, count, Especie)[,2] > 10

rep.sp <- ldply(Sp.raw.main.data, DoubleRep, detailed = TRUE)
rep.gen <- ldply(Gen.raw.main.data, DoubleRep, detailed = FALSE) 
rep.sp <- rep.sp[,-2]
REP <- rep.sp

names(REP)[1]<- "Specie"
names(REP)[2]<- "Trait"
names(REP)[3]<- "Repetability"
REP$Diff <- abs(REP$Repetability - REP$all.rep)

names(REP)

REP[,2] <- c( names(raw.main.data[[1]][12:50]) )
REP[,2] <- factor (REP[,2], levels = unique(REP[,2]) )

####### Plotando os valores de repetibilidade

REP<- tbl_df(REP)
REP %>% na.omit(.) %>%
  ggplot(., aes_string( x= "Trait", y="Repetability", color="Specie")) +
  geom_line(aes(group = Specie), size=2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

REP %>% na.omit(.) %>%
  ggplot(., aes_string( x= "Trait", y="all.rep" )) +
  #geom_line(aes(group = Specie)) +
  geom_point(aes(group = Specie),  size=1 ) +
  facet_wrap(~Specie, scale="free_y", nrow = 5, ncol = 7) +
  theme_bw() +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, size = 4), 
        axis.text.y = element_text(size = 7),
        axis.title.x = element_blank(),
        strip.text= element_text(size=7)) +
  ggtitle("Traits Repetabilities by Specie") + 
  theme(plot.title = element_text(lineheight=.8, face="bold", size = 8)) 

REP %>% na.omit(.) %>%
  ggplot(., aes_string( x= "Trait", y="Repetability" )) +
  #geom_line(aes(group = Specie)) +
  geom_point(aes(group = Specie),  size=1 ) +
  facet_wrap(~Specie, scale="free_y", nrow = 5, ncol = 7) +
  theme_bw() +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, size = 4), 
        axis.text.y = element_text(size = 7),
        axis.title.x = element_blank(),
        strip.text= element_text(size=7)) +
  ggtitle("Traits Repetabilities by Specie") + 
  theme(plot.title = element_text(lineheight=.8, face="bold", size = 8)) 


REP %>% na.omit(.) %>%
  ggplot(., aes_string( x= "Trait", y="Diff" )) +
  #geom_line(aes(group = Specie)) +
  geom_point(aes(group = Specie),  size=1 ) +
  facet_wrap(~Specie, scale="free_y", nrow = 5, ncol = 7) +
  theme_bw() +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, size = 4), 
        axis.text.y = element_text(size = 7),
        axis.title.x = element_blank(),
        strip.text= element_text(size=7)) +
  ggtitle("Traits Repetabilities by Specie") + 
  theme(plot.title = element_text(lineheight=.8, face="bold", size = 8)) 

####### wraped por caractere e com o nome das especies
REP %>% na.omit(.) %>%
  ggplot(., aes_string(x=1, y="Repetability" )) +
  #geom_line(aes(group = Specie)) +
  geom_point(aes(group = "Trait"),  size=1 ) +
  geom_text(aes(group = "Trait", label = Specie)) +
  facet_wrap(~Trait, scale="fixed", ncol =8, nrow = 5) +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=17),
        strip.text= element_text(size=15),
        plot.title = element_text(lineheight=.8, face="bold", size = 8)) 

############# de densidade
REP %>% na.omit(.) %>%
  ggplot(., aes_string(x = "Repetability" )) +
  geom_density(aes(group = Trait)) +
  #geom_text(aes(group = "Trait", label = Specie)) +
  #facet_wrap(~Trait, scale="free_y", ncol =5, nrow = 8) +
  theme_bw() +
  ylab("Density") +
  ggtitle("Density of traits' Repetabilities") + 
theme(legend.position="none",
      axis.text.x = element_blank(), 
      axis.text.y = element_text(size = 16),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=19),
      strip.text= element_text(size=15),
      plot.title = element_text(lineheight=.8, face="bold", size = 17))  

table(na.omit(REP$Repetability) >=0.95)

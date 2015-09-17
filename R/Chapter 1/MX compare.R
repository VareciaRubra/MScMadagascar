

  current.data <- sp.main.data
  current.data <- gen.main.data
  
  cov.mx <- current.data %>% llply(function(x) x$matrix$cov)
  cor.mx <- current.data %>% llply(function(x) x$matrix$cor)
  mx.rep <- current.data %>% ldply(function(x) x$Mx.Rep$BootsRep) 
  n.size <- current.data %>% ldply(function(x) x$sample.size) 
  
  boot.R2 <- current.data %>% llply(function(x) x$BootsR2)
  # criando mascaras para selecionar só parte do dataset.
  #Todas que tem matriz, mesmo as mal estimadas
  mask <- ldply(cov.mx, function(x) !is.na(x[1]))[,2]
  #Só mtrizes com mais de 40 indivíduos
  mask.n.size <- n.size[,2]>40
  
  Saguinus_G.cov <- read.csv(file = "Data/Saguinus_G.csv", header = F)
  Saguinus_G.cov <- as.matrix(Saguinus_G.cov)
  dimnames(Saguinus_G.cov) <- dimnames(cov.list[[1]])
  Saguinus_G.cor <- cov2cor(Saguinus_G.cov)
  
  Saguinus_P.cov <- read.csv(file = "Data/Saguinus_P.csv", header = F)
  Saguinus_P.cov <- as.matrix(Saguinus_P.cov)
  Saguinus_P.cov.t <- t(Saguinus_P.cov)
  Saguinus_P.cov[upper.tri(diag(39))] <- Saguinus_P.cov.t[upper.tri(diag(39))]
  
  dimnames(Saguinus_P.cov) <- dimnames(cov.list[[1]])
  Saguinus_P.cor <- cov2cor(Saguinus_P.cov)
  
  cov.list <- cov.mx[mask]
  cov.list$Saguinus_P.cov <- Saguinus_P.cov
  cov.list$Saguinus_G.cov <- Saguinus_G.cov
  
  cor.list <- cor.mx[mask]
  cor.list$Saguinus_P.cor <- Saguinus_P.cor
  cor.list$Saguinus_G.cor <- Saguinus_G.cor
  
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

  sample.size.list <- c(n.size[mask,2], 130, 230)
    
  
  mx.compare = vector("list", 5)
  mx.compare[1:5] <- NA
  mx.compare[[1]] <- RandomSkewers(cov.x= cov.list, num.vectors = 1000, repeat.vector = rep.list[, 2] )
  mx.compare[[2]]$correlations <- as.matrix(KrzCor(cov.x= cov.list, num.vectors = 1000, repeat.vector = rep.list[, 3] ))
  mx.compare[[3]]$correlations <- as.matrix(PCAsimilarity(cov.x= cov.list, num.vectors = 1000, repeat.vector = rep.list[, 4] ) )
  mx.compare[[4]]$correlations <- as.matrix(MatrixCor(cor.x= cor.list, correlation = TRUE, num.vectors = 1000, repeat.vector = rep.list[, 5] ) )
  mx.compare[[5]]$correlations <- as.matrix(KrzCor(cov.x= cor.list, correlation = TRUE, num.vectors = 1000, repeat.vector = rep.list[, 6] ))
  names(mx.compare)[1:5] <-  c('RS', 'KRZ','PCA.s', 'Mantel', 'KRZ')
  mx.class<- c('V/CV', 'V/CV','V/CV', 'COR', 'COR')
  for (i in 1:5)  {mx.compare[[i]]$method <- names(mx.compare)[i]}
  for (i in 1:5)  {mx.compare[[i]]$mx.class <- mx.class[i]}
  
  
  mat_data <- mx.compare$RS$correlations
  mat_data[lower.tri(mat_data)] <- t(mx.compare$KRZ$correlations)[lower.tri(mx.compare$KRZ$correlations)]
  #diag(mat_data) <- NA
  
  
  
   m.rs.krz = melt(mat_data) 
   m.rs.krz$Var1<- factor( m.rs.krz$Var1, levels = levels( m.rs.krz$Var1)[44:1])
   m.rs.krz.position =  m.rs.krz
   m.rs.krz.position$Var1 <- as.numeric( m.rs.krz.position$Var1)
   m.rs.krz.position$Var2 <- as.numeric( m.rs.krz.position$Var2)
   m.rs.krz.position$value= round( m.rs.krz.position$value, 3)
   
   myPalette <- colorRampPalette(brewer.pal(11, 'Spectral'), space = 'Lab')(n = 100)
   matrix_comparisons <- ggplot ( m.rs.krz) +
    geom_tile(aes(x = Var2, y = Var1, fill = value)) +
    scale_fill_gradientn(name = '', colours = myPalette) +
    ylab ('') + xlab ('') + labs(title = "Matrix comparisons") + 
    geom_text(data =  m.rs.krz.position, size = 3, aes(x = Var2, y = Var1, label = round(value, 2)) ) + 
    theme(axis.text.x = element_blank(),
          #axis.text.y = element_blank(),
          axis.ticks = element_line(size = 0),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 7),
          rect = element_blank(), line = element_blank())
  
  log.cov.mx <- current.data %>% llply(function(x) x$matrix$cov.log)
  log.cor.mx <- llply(log.cov.mx[mask], cov2cor)
  mx.compare.log = vector("list", 5)
  mx.compare.log[1:5] <- NA
  mx.compare.log[[1]] <- RandomSkewers(cov.x= log.cov.mx[mask], num.vectors = 1000, repeat.vector = mx.rep[mask, 2] )
  mx.compare.log[[2]]$correlations <- as.matrix(KrzCor(cov.x= log.cov.mx[mask], num.vectors = 1000, repeat.vector = mx.rep[mask, 3] ))
  mx.compare.log[[3]]$correlations <- as.matrix(PCAsimilarity(cov.x= log.cov.mx[mask], num.vectors = 1000, repeat.vector = mx.rep[mask, 4] ) )
  mx.compare.log[[4]]$correlations <- as.matrix(MatrixCor(cor.x= log.cor.mx, correlation = TRUE, num.vectors = 1000, repeat.vector = mx.rep[mask, 5] ) )
  mx.compare.log[[5]]$correlations <- as.matrix(KrzCor(cov.x= log.cor.mx, correlation = TRUE, num.vectors = 1000, repeat.vector = mx.rep[mask, 6] ))
  names(mx.compare.log)[1:5] <-  c('RS', 'KRZ','PCA.s', 'Mantel', 'KRZ')
  mx.class<- c('V/CV', 'V/CV','V/CV', 'COR', 'COR')
  for (i in 1:5)  {mx.compare.log[[i]]$method <- names(mx.compare.log)[i]}
  for (i in 1:5)  {mx.compare.log[[i]]$mx.class <- mx.class[i]}
  
  
##################################
#Função para plotar matrizes 
#com respectivos valores seguindo 
# palleta de padrao de cor
###################################
myPalette <- colorRampPalette(c("yellow", "white", "purple"))(n = 100)
myPalette <- "RdBu"  
plot.matrix<- function(mx = NULL, sample.size = NULL, brewer = myPalette, show.values = TRUE)
{ require (plotrix)
  library("RColorBrewer")
  paleta.of.choice  <- rev(brewer.pal(11, brewer))
  colfunc <- colorRampPalette(paleta.of.choice)
  paleta  <- colfunc(90)
  nCores  <- length(paleta)
  intervalMarks  <- seq(from = 0.1, to = 1, length.out = nCores)
  
  mx.cor <- mx
  mx.dimentions <- dim(mx.cor)[1]
  
  dados<- as.vector(mx.cor)
  intervals  <- findInterval(dados , intervalMarks, rightmost.closed = T, all.inside = T)
  cores  <- paleta[intervals]
  
  color2D.matplot(x = mx.cor, axes = F, cellcol = cores, show.values= 2, vcex= 0.6, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  #title(main = paste(mx$mx.class, "Matrices compared by", mx$method, "method"))
  title(main = "Matrix comparison via KRZ (lower) and RS (upper)")
  axis(1, 1:mx.dimentions, paste("(n=", sample.size,")",1:length(rownames(mx.cor)), sep = " "), las = 2, cex.axis = 0.9, tick = FALSE, line = 0)
  #axis(2, mx.dimentions:1, length(rownames(mx.cor)):1, las = 1, cex.axis = 1, tick = FALSE, line = 0)
  axis(2, mx.dimentions:1, paste(rownames(mx.cor), 1:length(rownames(mx.cor))), las = 1, cex.axis = 0.9, tick = FALSE, line = 0)

}
par(mar = c(7,17,3,2))
par(mfrow = c(3,1))
  lapply(mx.compare, plot.matrix)
  lapply(mx.compare.log, plot.matrix)
  
  
  plot.matrix(mx = mat_data, sample.size = sample.size.list, show.values = TRUE)
  #Calculando a similaridade média das matrizes.
mean.sim<-function(x) {
    x[upper.tri(x)]<-t(x)[upper.tri(x)]
    diag(x)<-NA
    x<-rowMeans(x,na.rm=TRUE)
    return(x)
}

##########Plotando os valores médios de comparaçao por metodo por matriz
current.mx.list <- mx.compare
current.mx.list <- mx.compare.log

mean.comp.values.cov<- current.mx.list[c(1:3)] %>% ldply (function (x) mean.sim(x$correlations) )
mean.comp.values.cor<- current.mx.list[c(4:5)] %>% ldply (function (x) mean.sim(x$correlations) )
mean.comp.values <- mean.comp.values.cov
names(mean.comp.values)[1] <- "method"

method.mean <- data.frame("method"  = mean.comp.values[,1], "mean"=apply(mean.comp.values[,-1], 1, mean))
Mean.Mx.Plots <- mean.comp.values %>% 
  gather(key=.sp, value = value, 2:43) %>%
  ggplot( .,aes(x= method, y = value, shape = method), varwidth = T) +
  geom_violin(aes(label = .sp, color = method, shape = method), alpha = 0.5) +
  geom_text( aes(label = .sp), size =2, vjust = 1, alpha = 0.4)  +
  geom_jitter(aes(shape = method, color = method)) +
  #geom_line(aes(group = .sp, color =.sp)) +
  scale_shape( guide = "none", name = "Mean value by method") +
  ggtitle("Mean value of comparisson") +
  theme_bw() +
  theme(plot.title = element_text(lineheight=.8, face="bold"))

Mean.Mx.Plots + geom_point(aes (x=method, y=mean, color = method, shape = method), data = method.mean, size = 6 )


########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
PCs1to4<- current.data[mask] %>% llply(function(x) as.data.frame(eigen(x$matrix$cov)$vectors[,1:4]) ) %>% llply(function(x) as.list(x) )
PC1 <- current.data[mask] %>% ldply(function(x) eigen(x$matrix$cov)$vectors[,1])
names.sp <-  PC1[,1]
PC1<- as.matrix(PC1[,-1])
dimnames(PC1)[[1]] <- names.sp
# se estiver correlacionado negativamente com o vetor isometrico, multiplicar por -1
ed.names <- names(current.data[mask][[1]]$ed)
Iso.Compare <- function(x) {
  isometrico<- rep( (1/sqrt(39)) , 39)
  corr = rep(NA, length(x) )
  re.oriented = x
    for (i in 1:length(x)) {
    corr[i]<- x[[i]] %*% isometrico
      
      if(corr[i] <0) {re.oriented[[i]] <- as.numeric(as.matrix(x[[i]]) %*% -1)}  else x[[i]]
    }
    
  names (corr) <- paste ("PC", 1:length(x), sep = "")
  names (re.oriented) <- paste ("PC", 1:length(x), sep = "") 
  
  return(results= list("corr" = corr, "re.oriented" = re.oriented)) }
PCs1to4<- current.data[mask] %>% llply(function(x) as.data.frame(eigen(x$matrix$cov)$vectors[,1:4]) ) %>% llply(function(x) as.list(x) )
PCs1to4.log<- current.data[mask] %>% llply(function(x) as.data.frame(eigen(x$matrix$cov.log)$vectors[,1:4]) ) %>% llply(function(x) as.list(x) )

Iso.Compare.reoriented <- llply(PCs1to4.log, Iso.Compare) %>% ldply(function(x) as.data.frame(x$re.oriented)) 

names(Iso.Compare.reoriented) <- c(".sp", "PC1", "PC2", "PC3", "PC4")
Iso.Compare.reoriented$.ed <- ed.names
Iso.Compare.reoriented %>% gather(key=".pcScore", value=value, 2:5 )
Iso.Compare.reoriented$.sp <- factor(Iso.Compare.reoriented$.sp, levels = unique(Iso.Compare.reoriented$.sp)) 
Iso.Compare.reoriented$.ed <- factor(Iso.Compare.reoriented$.ed, levels = rev(ed.names) ) 
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
Iso.Compare.reoriented %>% gather(key=".pcScore", value=value, 2:5 ) %>%
  ggplot (.) +
  geom_tile(aes(x = .pcScore, y = .ed, fill = value)) +
  facet_wrap(~.sp) +
  theme_bw() +
  scale_fill_gradientn(name = 'PC.Score', colours = myPalette(1000)) +
  ylab ('PC.Scores') + xlab ('') + labs(title = "First 4 PC's scores") +
  #scale_y_discrete(limits = levels(PCs1to4$.pcScore))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks = element_line(size = 0),
        axis.text.y = element_text(size= 7),
        strip.text= element_text(size=7)) +
  coord_flip()

Iso.Compare.cor.iso <- llply(PCs1to4, Iso.Compare) %>% ldply(function(x) x$corr) 
names(Iso.Compare.cor.iso) <- c(".sp", "PC1", "PC2", "PC3", "PC4")
Iso.Compare.cor.iso$.sp <- factor(Iso.Compare.cor.iso$.sp, levels = rev(unique(Iso.Compare.cor.iso$.sp))) 
Iso.Compare.cor.iso %>% gather(key="Isometric.Correlation", value=value, 2:5 ) %>%
  ggplot (.) +
  geom_tile(aes(x = Isometric.Correlation, y = .sp, fill = abs(value) ) ) +
  theme_bw() +
  scale_fill_gradientn(name = 'Isometric.Correlation', colours = myPalette(1000)) +
  ylab ('') + xlab ('') + labs(title = "Matrix of log values correlation with isometric vector") +
  #scale_y_discrete(limits = levels(PCs1to4$.pcScore))) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10),
        axis.text.x = element_text(size = 10),
        axis.ticks = element_line(size = 0),
        axis.text.y = element_text(size= 10)) +
  theme(plot.title = element_text(lineheight=.8, face="bold"))

PCs1to4<- current.data[mask] %>% ldply(function(x) as.data.frame(eigen(x$matrix$cov)$vectors[,1:4]) )  
names(PCs1to4) <- c(".sp", "PC1", "PC2", "PC3", "PC4")
PCs1to4$.ed <- ed.names
PCs1to4 %>% gather(key=".pcScore", value=value, 2:5 )
PCs1to4$.sp <- factor(PCs1to4$.sp, levels = unique(PCs1to4$.sp)) 
PCs1to4$.ed <- factor(PCs1to4$.ed, levels = rev(ed.names) ) 
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
PCs1to4 %>% gather(key=".pcScore", value=value, 2:5 ) %>%
  ggplot (.) +
  geom_tile(aes(x = .sp, y = .ed, fill = value)) +
  facet_wrap(~.pcScore) +
  theme_bw() +
  scale_fill_gradientn(name = 'PC.Score', colours = myPalette(1000)) +
  ylab ('PC.Scores') + xlab ('') + labs(title = "First 4 PC's scores") +
  #scale_y_discrete(limits = levels(PCs1to4$.pcScore))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks = element_line(size = 0))

data.frame ('n.size' = n.size [mask, -1], Iso.Compare.cor.iso) %>%
  gather(key = '.pc', value=value, c(3:6)) %>%
  ggplot (.) +
  #geom_text(aes (x = n.size, y = abs(value), label = .sp)) +
  geom_point(aes (x = n.size, y = abs(value), size = n.size, color = abs(value) ) ) +
  geom_text(aes (x = n.size, y = abs(value), label = .sp), alpha = 0.4, size =2) +
  #scale_color_brewer(name = 'PC.Score' , type = "div"  , palette =  myPalette(1000)) +
  ylab ('Absolute value of correlation') + xlab ('Sample size') + labs(title = "First 4 PC's correlation with Isometric vector sample size") +
  facet_wrap(~ .pc) +
  theme_bw() +
  theme(plot.title = element_text(lineheight=.8, face="bold"))

rarefaction.all.rs <- llply(sp.main.data, function (x) x$rarefacation$rs)

rarefaction.1 <-PlotRarefaction(sp.main.data$Euoticus_elegantulus$rarefaction$rs) + labs(title = "Euoticus senegalensis") 
rarefaction.2 <-PlotRarefaction(sp.main.data$Loris_tardigradus$rarefaction$rs) + labs(title = "Loris tardigradus") 
rarefaction.3 <-PlotRarefaction(sp.main.data$Daubentonia_madagascariensis$rarefaction$rs) + labs(title = "D. madagascariensis") 
rarefaction.4 <-PlotRarefaction(sp.main.data$Phaner_furcifer$rarefaction$rs) + labs(title = "Phaner furcifer") 
rarefaction.5 <-PlotRarefaction(sp.main.data$Microcebus_griseorufus$rarefaction$rs) + labs(title = "Microcebus griseorufus") 
rarefaction.6 <-PlotRarefaction(sp.main.data$Indri_indri$rarefaction$rs) + labs(title = "Indri indri") 



multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


multiplot(rarefaction.5, rarefaction.6, cols = 2)


########################################################################################################################################
########################################################################################################################################
########################################################################################################################################


########### quantidade de variação espalhada nos primeiros 4 PC
  PcPercent.1to4<- cov.mx[mask] %>% ldply(function(x) eigen(x)$values[1:4]/sum(eigen(x)$values) )
  names(PcPercent.1to4) <- c(".sp", "PC1", "PC2", "PC3", "PC4")
  PcPercent.1to4$.sp <- factor(PcPercent.1to4$.sp, levels = unique(PcPercent.1to4$.sp)) 
  library(scales)
  PcPercent.1to4 %>%
  gather(key="Percent_var_PC", value=value, 2:5 ) %>%
    ggplot( ., aes(x= Percent_var_PC, y = value, color = Percent_var_PC, label = .sp), varwidth = T) +
    scale_y_continuous(labels=percent) +
    geom_text(size = 2, vjust = 1, alpha = 0.4 )  +
    theme_bw() +
    geom_jitter() +
    geom_violin(alpha = 0) +
    ggtitle("% of variance on first 4 PC by specie") +
    theme(plot.title = element_text(lineheight=.8, face="bold"), axis.title.x = element_blank()) +
    theme(legend.position="none")

PcPercent.1to4<- cov.mx[mask] %>% ldply(function(x) eigen(x)$values[1:4]/sum(eigen(x)$values) )
names(PcPercent.1to4) <- c(".sp", "PC1", "PC2", "PC3", "PC4")
PcPercent.1to4$.sp <- factor(PcPercent.1to4$.sp, levels = unique(PcPercent.1to4$.sp)) 
library(scales)
PcPercent.1to4 %>%
  gather(key="Percent_var_PC", value=value, 2:5 ) %>%
  ggplot( ., aes(x= Percent_var_PC, y = value, color = .sp), varwidth = T) +
  geom_line(aes(group = .sp), size = 1) +
  scale_y_continuous(labels=percent) +
  theme_bw() +
  facet_wrap(~.sp)+
  ggtitle("% of variance on first 4 PC by specie") +
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  theme(axis.title.x = element_blank()) +
  theme(legend.position="none")

MMxStats<- cov.mx[mask] %>% ldply(function(x) MeanMatrixStatistics(x))
names(MMxStats)[1] <- ".sp"
names(MMxStats)[3] <- "PC1.percent"
MMxStats %>% 
gather(key= .MMxStats, value=value, 2:10 ) %>%
ggplot( ., aes(x= .MMxStats, y = value, color = .MMxStats, label = .sp), varwidth = T) +
  geom_text(size =2, vjust = 1)  +
  theme_bw() +
  geom_boxplot(alpha = 0) +
  facet_wrap(~.MMxStats, scales="free") +
  ggtitle("Mean matrix evolutionary statistics by specie") +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank()) +
  theme(plot.title = element_text(lineheight=.8, face="bold")) + 
  theme(legend.position="none")


MMxStats<- cov.list %>% ldply(function(x) MeanMatrixStatistics(x))
names(MMxStats)[1] <- ".sp"
MMxStats[,1] <- factor(unique(MMxStats[,1]), levels = unique(MMxStats[,1]))
names(MMxStats)[2] <- "R^2"
names(MMxStats)[3] <- "PC1.percent"

MMxStats %>% .[, c(1:3,9)] %>%
  gather(key= .MMxStats, value=value, 2:4 ) %>% 
  ggplot(aes(y=  value, x = .sp, fill = .MMxStats ), varwidth = T) +
  geom_bar(aes(y=  value, x = .sp, fill = .MMxStats ), stat = "identity", position="dodge")+
  #scale_y_continuous(breaks = seq(0,0.6, by = 0.1) )+
  theme_bw() +
  facet_grid(~.MMxStats) +
  coord_flip() +
  ggtitle("Mean matrix evolutionary statistics by specie") +
    theme(axis.text.y = element_text(face = "italic")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  scale_fill_grey() +
  #scale_fill_brewer(palette="Spectral") +
  scale_colour_brewer(name = "Mean Matrix Statistics") + 
  theme(legend.position="none")




MMxStats<- cov.mx[mask] %>% ldply(function(x) MeanMatrixStatistics(x))
names(MMxStats)[1] <- ".sp"
MMxStats %>% 
  ggplot( ., aes(x = MeanSquaredCorrelation, y = ICV, color = .sp, label = .sp), varwidth = T) +
  geom_text(size =4, vjust = 1 )  +
  stat_smooth(method="lm", aes(group=1)) +
  scale_x_continuous(limits = c(0.03, 0.4)) +
  scale_y_continuous(limits = c(1.5, 4.5)) +
  theme_bw() +
  geom_point() +
  #geom_text(x = 0.5, y = 3, label = "r.squared = 0.899 \n 38 DF,  p-value: < 2.2e-16") + 
  ggtitle("R2 x ICV") +
  theme(legend.position  ="none") +
  theme(plot.title = element_text(lineheight=.8, face="bold"))
summary(lm(ICV ~ MeanSquaredCorrelation, data = MMxStats))

MMxStats %>% 
  ggplot( ., aes(x = flexibility, y = evolvability, color = .sp, label = .sp), varwidth = T) +
  geom_text(size =5, vjust = 1) +
  stat_smooth(method="lm", aes(group=1)) +
  #scale_x_continuous(limits = c(0.23, 0.5)) +
  #scale_y_continuous(limits = c(1, 7)) +
  theme_bw() +
  geom_point() +
  #geom_text(x = 0.35, y = 3.5, label = "r.squared = 0.149 \n 38 DF,  p-value: 0.00803") + 
  ggtitle("flexibility x evolvability") +
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  theme(legend.position  ="none") 
summary(lm(flexibility ~ evolvability , data = MMxStats))

names (MMxStats)[3] <- "PC1.percent"
MMxStats %>% 
  ggplot( ., aes(x = flexibility, y = PC1.percent, color = .sp, label = .sp), varwidth = T) +
  geom_text(size =4, vjust = 1)  +
  stat_smooth(method="lm", aes(group=1)) +
  #scale_x_continuous(limits = c(0.1, 0.7)) +
  #scale_y_continuous(limits = c(1, 7)) +
  theme_bw() +
  geom_point() +
  #geom_text(x = 0.4,  y = 0.74, label = "r.squared = 0.8992 \n 38 DF,  p-value: < 2.2e-16") + 
  ggtitle("flexibility x PC1.percent") +
  theme(legend.position  ="none") +
  theme(plot.title = element_text(lineheight=.8, face="bold"))

summary(lm(flexibility ~ PC1.percent , data = MMxStats))

gm.mean <- current.data[mask] %>% ldply(function(x) x$gm.mean)

MMxStats.gm<- cbind(MMxStats, gm.mean[,2])
names(MMxStats.gm)[11] <- "gm.mean"

names (MMxStats.gm)[3] <- "PC1.percent"
MMxStats.gm %>%
  ggplot( ., aes(x = gm.mean, y = PC1.percent, color = .sp, label = .sp), varwidth = T) +
  geom_text(size =5, vjust = 1)  +
  stat_smooth(method="lm", aes(group=1)) +
  #scale_x_continuous(limits = c(0.1, 0.7)) +
  #scale_y_continuous(limits = c(1, 7)) +
  theme_bw() +
  geom_point() +
  #geom_text(x = 0.4,  y = 0.74, label = "r.squared = 0.8992 \n 38 DF,  p-value: < 2.2e-16") + 
  ggtitle("gm.mean x PC1.percent") +
  theme(legend.position  ="none") +
  theme(plot.title = element_text(lineheight=.8, face="bold"))

summary(lm(gm.mean ~ PC1.percent , data = MMxStats.gm))




  current.data <- sp.main.data
  current.data <- gen.main.data
  
  cov.mx <- current.data %>% llply(function(x) x$matrix$cov)
  cor.mx <- current.data %>% llply(function(x) x$matrix$cor)
  mx.rep <- current.data %>% ldply(function(x) x$Mx.Rep$BootsRep) 
  n.size <- current.data %>% ldply(function(x) x$sample.size) 
  # criando mascaras para selecionar só parte do dataset.
  #Todas que tem matriz, mesmo as mal estimadas
  mask <- ldply(cov.mx, function(x) !is.na(x[1]))[,2]
  #Só mtrizes com mais de 40 indivíduos
  mask.n.size <- n.size[,2]>40
  
  mx.compare = vector("list", 5)
  mx.compare[1:5] <- NA
  mx.compare[[1]] <- RandomSkewers(cov.x= cov.mx[mask], num.vectors = 1000, repeat.vector = mx.rep[mask, 2] )
  mx.compare[[2]]$correlations <- as.matrix(KrzCor(cov.x= cov.mx[mask], num.vectors = 1000, repeat.vector = mx.rep[mask, 3] ))
  mx.compare[[3]]$correlations <- as.matrix(PCAsimilarity(cov.x= cov.mx[mask], num.vectors = 1000, repeat.vector = mx.rep[mask, 4] ) )
  mx.compare[[4]]$correlations <- as.matrix(MatrixCor(cor.x= cor.mx[mask], correlation = TRUE, num.vectors = 1000, repeat.vector = mx.rep[mask, 5] ) )
  mx.compare[[5]]$correlations <- as.matrix(KrzCor(cov.x= cor.mx[mask], correlation = TRUE, num.vectors = 1000, repeat.vector = mx.rep[mask, 6] ))
  names(mx.compare)[1:5] <-  c('RS', 'KRZ','PCA.s', 'Mantel', 'KRZ')
  mx.class<- c('V/CV', 'V/CV','V/CV', 'COR', 'COR')
  for (i in 1:5)  {mx.compare[[i]]$method <- names(mx.compare)[i]}
  for (i in 1:5)  {mx.compare[[i]]$mx.class <- mx.class[i]}
  
  
##################################
#Função para plotar matrizes 
#com respectivos valores seguindo 
# palleta de padrao de cor
###################################
plot.matrix<- function(mx = NULL, brewer = "RdBu", show.values = TRUE)
{ require (plotrix)
  library("RColorBrewer")
  paleta.of.choice  <- rev(brewer.pal(11, brewer))
  colfunc <- colorRampPalette(paleta.of.choice)
  paleta  <- colfunc(90)
  nCores  <- length(paleta)
  intervalMarks  <- seq(from = 0.1, to = 1, length.out = nCores)
  
  mx.cor <- mx$correlations
  mx.dimentions <- dim(mx.cor)[1]
  
  dados<- as.vector(mx.cor)
  intervals  <- findInterval(dados , intervalMarks, rightmost.closed = T, all.inside = T)
  cores  <- paleta[intervals]
  
  color2D.matplot(x = mx.cor, axes = F, cellcol = cores, show.values= 2, vcex= 0.4, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  title(main = paste(mx$mx.class, "Matrices compared by", mx$method, "method"))
  axis(1, 1:mx.dimentions, paste("(n=", n.size[mask,2],")",1:length(rownames(mx.cor)), sep = " "), las = 2, cex.axis = 0.7, tick = FALSE, line = 0)
  #axis(2, mx.dimentions:1, length(rownames(mx.cor)):1, las = 1, cex.axis = 1, tick = FALSE, line = 0)
  axis(2, mx.dimentions:1, paste(rownames(mx.cor), 1:length(rownames(mx.cor))), las = 1, cex.axis = 0.7, tick = FALSE, line = 0)

}
par(mar = c(5,15,3,2))
par(mfrow = c(3,1))
  lapply(mx.compare, plot.matrix)
#Calculando a similaridade média das matrizes.
mean.sim<-function(x) {
    x[upper.tri(x)]<-t(x)[upper.tri(x)]
    diag(x)<-NA
    x<-rowMeans(x,na.rm=TRUE)
    return(x)
}

##########Plotando os valores médios de comparaçao por metodo por matriz
lapply(mx.compare, mean.sim)
mean.comp.values.cov<- mx.compare[c(1:3)] %>% ldply (function (x) mean.sim(x$correlations) )
mean.comp.values.cor<- mx.compare[c(4:5)] %>% ldply (function (x) mean.sim(x$correlations) )
mean.comp.values <- mean.comp.values.cor
names(mean.comp.values)[1] <- "method"

method.mean <- data.frame("method"  = mean.comp.values[,1], "mean"=apply(mean.comp.values[,-1], 1, mean))
Mean.Mx.Plots <- mean.comp.values %>% 
  gather(key=.sp, value = value, 2:43) %>%
  ggplot( .,aes(x= method, y = value, shape = method), varwidth = T) +
  geom_violin(aes(label = .sp, color = method, shape = method), alpha = 0.5) +
  geom_text( aes(label = .sp), size =5, vjust = 1, alpha = 0.4)  +
  geom_jitter(aes(shape = method, color = method)) +
  #geom_line(aes(group = .sp, color =.sp)) +
  scale_shape( guide = "none", name = "Mean value by method") +
  ggtitle("Mean value of comparisson for each species and by method") +
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
    
  #rownames (corr) <- paste ("PC", 1:length(x), sep = "")
  names (re.oriented) <- paste ("PC", 1:length(x), sep = "") 
  
  return(results= list("corr" = corr, "re.oriented" = re.oriented)) }
PCs1to4<- current.data[mask] %>% llply(function(x) as.data.frame(eigen(x$matrix$cov)$vectors[,1:4]) ) %>% llply(function(x) as.list(x) )
Iso.Compare.reoriented <- llply(PCs1to4, Iso.Compare) %>% ldply(function(x) as.data.frame(x$re.oriented)) 
names(Iso.Compare.reoriented) <- c(".sp", "PC1", "PC2", "PC3", "PC4")
Iso.Compare.reoriented$.ed <- ed.names
Iso.Compare.reoriented %>% gather(key=".pcScore", value=value, 2:5 )
Iso.Compare.reoriented$.sp <- factor(Iso.Compare.reoriented$.sp, levels = unique(Iso.Compare.reoriented$.sp)) 
Iso.Compare.reoriented$.ed <- factor(Iso.Compare.reoriented$.ed, levels = rev(ed.names) ) 
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
Iso.Compare.reoriented %>% gather(key=".pcScore", value=value, 2:5 ) %>%
  ggplot (.) +
  geom_tile(aes(x = .sp, y = .ed, fill = value)) +
  facet_wrap(~.pcScore) +
  theme_bw() +
  scale_fill_gradientn(name = 'PC.Score', colours = myPalette(1000)) +
  ylab ('PC.Scores') + xlab ('') + labs(title = "First 4 PC's scores") +
  #scale_y_discrete(limits = levels(PCs1to4$.pcScore))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks = element_line(size = 0),
        axis.text.y = element_text(size= 7)) 

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
    geom_text(size = 5, vjust = 1 )  +
    theme_bw() +
    geom_jitter() +
    geom_violin(alpha = 0) +
    ggtitle("% of variance on first 4 PC by specie") +
    theme(plot.title = element_text(lineheight=.8, face="bold"))

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
  geom_text(size =5, vjust = 1)  +
  theme_bw() +
  geom_boxplot(alpha = 0) +
  facet_wrap(~.MMxStats, scales="free") +
  ggtitle("Mean matrix evolutionary statistics by specie") +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank()) +
  theme(plot.title = element_text(lineheight=.8, face="bold")) + 
  theme(legend.position="none")

MMxStats<- cov.mx[mask] %>% ldply(function(x) MeanMatrixStatistics(x))
names(MMxStats)[1] <- ".sp"
MMxStats %>% 
  ggplot( ., aes(x = MeanSquaredCorrelation, y = ICV, color = .sp, label = .sp), varwidth = T) +
  geom_text(size =4, vjust = 1 )  +
  stat_smooth(method="lm", aes(group=1)) +
  scale_x_continuous(limits = c(0.03, 0.6)) +
  scale_y_continuous(limits = c(1.5, 5)) +
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
  scale_x_continuous(limits = c(0.23, 0.5)) +
  #scale_y_continuous(limits = c(1, 7)) +
  theme_bw() +
  geom_point() +
  #geom_text(x = 0.35, y = 3.5, label = "r.squared = 0.149 \n 38 DF,  p-value: 0.00803") + 
  ggtitle("flexibility x evolvability") +
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  theme(legend.position  ="none") 
summary(lm(flexibility ~ evolvability , data = MMxStats))

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

MMxStats.gm %>% 
  #filter(.sp != "Cheirogaleus_major") %>%
  #filter(.sp != "Cheirogaleus_medius") %>%
  #filter(.sp != "Eulemur_rufifrons") %>%
  #filter(.sp != "Phaner_furcifer") %>%
  #filter(.sp != "Loris_tardigradus") %>%
  #filter(.sp != "Lepilemur_leucopus") %>%
  #filter(.sp != "Loris_tardigradus") %>%
  #filter(PC1.percent <=0.4)%>%
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


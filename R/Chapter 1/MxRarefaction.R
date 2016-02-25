mean.sim<-function(x) {
    x[lower.tri(x)]<-t(x)[lower.tri(x)]
  diag(x)<-NA
  x<-rowMeans(x, na.rm=TRUE)
  return(x)
}

min.sim<-function(x) {
  x[lower.tri(x)]<-t(x)[lower.tri(x)]
  diag(x)<-NA
  t_<- dim(x)[[1]]
  x<- apply(x, 2, FUN = min, na.rm = T )
  return(x)
}

max.sim<-function(x) {
  x[lower.tri(x)]<-t(x)[lower.tri(x)]
  diag(x)<-NA
  x<-apply(x, 2, FUN = max, na.rm = T)
  return(x)
}

ic.sim<-function(x) {
  x[lower.tri(x)]<-t(x)[lower.tri(x)]
  diag(x)<-NA
  t_<- dim(x)[[1]]
  ic<- apply(x, 2, FUN = function (x) sort(x)[round(t_*c(0.025, 0.975) )]  )
  names(ic) <- c("min", "max")
  
  return(ic)
}

ic.max.sim<-function(x) {
  x[lower.tri(x)]<-t(x)[lower.tri(x)]
  diag(x)<-NA
  t_<- dim(x)[[1]] -1
  ic<- apply(x, 2, FUN = function (a) sort(a)[t_-1]  )
  return(ic)
}

ic.min.sim<-function(x) {
  x[lower.tri(x)]<-t(x)[lower.tri(x)]
  diag(x)<-NA
  t_<- dim(x)[[1]]
  ic<- apply(x, 2, FUN = function (a) sort(a)[2]  )
  return(ic)
}

mean.compare.cov <- mx.compare[c(1:2)] %>% llply (function (x) x$correlations[-c(43:44), -c(43:44)] ) %>% ldply (function (x) mean.sim(x) )
min.compare.cov <- mx.compare[c(1:2)] %>% llply (function (x) x$correlations[-c(43:44), -c(43:44)] ) %>% ldply (function (x) min.sim(x) )
max.compare.cov <- mx.compare[c(1:2)] %>% llply (function (x) x$correlations[-c(43:44), -c(43:44)] ) %>% ldply (function (x) max.sim(x) )
ic.compare.cov <- mx.compare[c(1:2)] %>% llply (function (x) x$correlations[-c(43:44), -c(43:44)] ) %>% ldply (function (x) ic.sim(x) )
ic.min.compare.cov <- mx.compare[c(1:2)] %>% llply (function (x) x$correlations[-c(43:44), -c(43:44)] ) %>% ldply (function (x) ic.min.sim(x) )
ic.max.compare.cov <- mx.compare[c(1:2)] %>% llply (function (x) x$correlations[-c(43:44), -c(43:44)] ) %>% ldply (function (x) ic.max.sim(x) )


dimnames(mean.compare.cov)[[1]] <- c(paste("Mean", mean.compare.cov[,1], sep = "."))
dimnames(min.compare.cov)[[1]] <- c(paste("Min", min.compare.cov[,1], sep = "."))
dimnames(max.compare.cov)[[1]] <- c(paste("Max", max.compare.cov[,1], sep = "."))
dimnames(ic.compare.cov)[[1]][c(1,3)] <- c(paste("ic.min", ic.compare.cov[,1][c(1,3)], sep = "."))
dimnames(ic.compare.cov)[[1]][c(2,4)] <- c(paste("ic.max", ic.compare.cov[,1][c(2,4)], sep = "."))
dimnames(ic.max.compare.cov)[[1]] <- c(paste("Max", ic.max.compare.cov[,1], sep = "."))
dimnames(ic.min.compare.cov)[[1]] <- c(paste("Min", ic.min.compare.cov[,1], sep = "."))



range.compare <- rbind(mean.compare.cov, ic.min.compare.cov, ic.max.compare.cov, ic.compare.cov)
range.compare <- range.compare[,-1]
range.compare <- as.data.frame(t(range.compare))

sample.size.list <- n.size[mask,2]
range.compare$n <- sample.size.list
range.compare$.id <- rownames(range.compare)

############# Rarefação no calculo das matrizes 
# pegando as maiores amostras de duas espécies nos extremos da variaçao de tamanho : Microcebus griseorrufus e Indri indri
# por dois métodos de comparação RS e KRZ 
# tudo ta calculado já no sp.main.data

PlotRaRangefaction<- function (comparison.list, y.axis = "Random Skewers", x.axis = "Number of sampled specimens") 
 ## comparison.list ==> output da evolqg::RarefactionStat()
  y.axis = "Random Skewers"
  x.axis = "Number of sampled specimens"
  comparison.list <- sp.main.data$Microcebus_griseorufus$rarefaction$rs
  plot.df <- melt(comparison.list, value.name = "Mean.BS.RS")
  names(plot.df) <- c("Mean.BS.RS", "n") 
  plot.df <- as.data.frame(lapply(plot.df, as.numeric))
rarefaction.plot.rs <- 
    ggplot(data = plot.df, aes(x= n, y = Mean.BS.RS)) + 
    geom_boxplot(aes( group = n), color="grey") + 
    scale_x_continuous(x.axis) + 
    geom_point(data = range.compare, aes(x = n, y = Mean.BS.RS), shape=17, size = 2.2, color = "red") +
    geom_text(data = range.compare, aes(x = n, y = Mean.BS.RS, label = .id), color = "red", size = 3, alpha = 0.3, angle = 70, hjust =  0.6) +
    geom_errorbar(data = range.compare, aes(x= n, ymin = Min.BS.RS, ymax = Max.BS.RS), alpha = 0.3, size = 1.2, color = "red") +
    geom_hline(data = range.compare, aes(yintercept = mean (Mean.BS.RS), y = mean (Mean.BS.RS)), alpha = 0.1, size = 2.5) +
    ggtitle("Rarefaction and Observed range of similarity") +
    theme(plot.title = element_text(lineheight=.8, face="bold"),
            axis.text.x = element_text(size = 15), 
            axis.title.x = element_text(size = 15)) +
    scale_y_continuous(y.axis, limits = c(0,1), breaks = c(0.00, 0.25, 0.5, 0.66, 0.75, 1.00)) +
    scale_x_continuous(x.axis, limits = c(0, 68)) +
    theme_bw()
    
    y.axis = "Krzanowski "
    x.axis = "Number of sampled specimens"
    comparison.list <- sp.main.data$Microcebus_griseorufus$rarefaction$krz
    plot.df <- melt(comparison.list, value.name = "Mean.BS.KRZ")
    names(plot.df) <- c("Mean.BS.KRZ", "n") 
    plot.df <- as.data.frame(lapply(plot.df, as.numeric))
rarefaction.plot.krz <- 
    ggplot(data = plot.df, aes(x= n, y = Mean.BS.KRZ)) + 
      geom_boxplot(aes( group = n), color="grey") + 
      scale_x_continuous(x.axis, name = "Number of sampled specimens") + 
      geom_point(data = range.compare, aes(x = n, y = Mean.BS.KRZ), shape=17, size = 2.2, color = "red") +
      geom_text(data = range.compare, aes(x = n, y = Mean.BS.KRZ, label = .id), color = "red", size = 3, alpha = 0.3, angle = 70, hjust =  0.6) +
      geom_errorbar(data = range.compare, aes(x= n, ymin = Min.BS.KRZ, ymax = Max.BS.KRZ), alpha = 0.3, size = 1.2, color = "red") +
      geom_hline(data = range.compare, aes(yintercept = mean (Mean.BS.KRZ), y = mean (Mean.BS.KRZ)), alpha = 0.1, size = 2.5) +
      ggtitle("Rarefaction and Observed range of similarity") +
      theme(plot.title = element_text(lineheight=.8, face="bold"),
            axis.text.x = element_text(size = 15), 
           labs(x = x.axis),
            axis.title.x = element_text(size = 15)) +
      scale_y_continuous(y.axis, limits = c(0,1), breaks = c(0.00, 0.25, 0.5, 0.75, 0.82, 1.00)) +
      scale_x_continuous(x.axis, limits = c(0, 68)) +
      theme_bw()


micro.rar.obs.sim <- plot_grid(rarefaction.plot.rs, rarefaction.plot.krz,labels = LETTERS[1:2])





### Add os ranges e valores médios de comparação por espécie no grafico de rarefaçao

p1 <- PlotRarefaction(sp.main.data$Microcebus_griseorufus$rarefaction$krz) + 
  labs(title = "Microcebus griseorufus") + 
  theme(plot.title = element_text(face = "italic", size = 20), 
        axis.title.x = element_blank(), 
        axis.text.x =element_blank(), 
        axis.ticks.x = element_line(size =0), 
        axis.text= element_text(size = 15)) + 
  geom_errorbar(data = range.compare, aes(x= L1, ymin = Min.BS.KRZ, ymax = Max.BS.KRZ)) + 
  geom_point(data = range.compare, aes(x = L1, y = Mean.BS.KRZ)) +
  coord_cartesian(ylim=c(0.4, 1)) + 
  scale_y_continuous("KRZ", breaks=seq(0.4, 1, 0.2))


p2 <- PlotRarefaction(sp.main.data$Indri_indri$rarefaction$krz) + 
  labs(title = "Indri indri") + 
  theme(plot.title = element_text(face = "italic", size = 20), 
        axis.title = element_blank(), 
        axis.text =element_blank(), 
        axis.ticks = element_line(size =0)) +
  coord_cartesian(ylim=c(0.4, 1)) + 
  scale_y_continuous( breaks=seq(0.40, 1, 0.2))

p3 <- PlotRarefaction(sp.main.data$Microcebus_griseorufus$rarefaction$rs) + 
  theme(axis.text= element_text(size = 15)) +
  coord_cartesian(ylim=c(0, 1)) + 
  scale_y_continuous("RS", breaks=seq(0.25, 1, 0.25))

p4 <- PlotRarefaction(sp.main.data$Indri_indri$rarefaction$rs) + 
  theme(axis.text= element_text(size = 15), 
        axis.title.y = element_blank(), 
        axis.text.y =element_blank(), 
        axis.ticks.y = element_line(size =0)) +
  coord_cartesian(ylim=c(0, 1)) + 
  scale_y_continuous(breaks=seq(0.25, 1, 0.25))

plot_grid(p1, p2, p3, p4)  






  
  



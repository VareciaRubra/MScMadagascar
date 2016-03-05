############# Rarefação no calculo das matrizes 
# pegando as maiores amostras de duas espécies nos extremos da variaçao de tamanho : Microcebus griseorrufus e Indri indri
# por dois métodos de comparação RS e KRZ 
# tudo ta calculado já no sp.main.data

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

### Add os ranges e valores médios de comparação por espécie no grafico de rarefaçao

mean.compare.cov <- mx.compare[c(1:2)] %>% llply (function (x) x$correlations[c(3:30,32:34), c(3:30,32:34)] ) %>% ldply (function (x) mean.sim(x) )
min.compare.cov <- mx.compare[c(1:2)] %>% llply (function (x) x$correlations[c(3:30,32:34), c(3:30,32:34)] ) %>% ldply (function (x) min.sim(x) )
max.compare.cov <- mx.compare[c(1:2)] %>% llply (function (x) x$correlations[c(3:30,32:34), c(3:30,32:34)] ) %>% ldply (function (x) max.sim(x) )
ic.compare.cov <- mx.compare[c(1:2)] %>% llply (function (x) x$correlations[c(3:30,32:34), c(3:30,32:34)] ) %>% ldply (function (x) ic.sim(x) )
ic.min.compare.cov <- mx.compare[c(1:2)] %>% llply (function (x) x$correlations[c(3:30,32:34), c(3:30,32:34)] )  %>% ldply (function (x) ic.min.sim(x) )
ic.max.compare.cov <- mx.compare[c(1:2)] %>% llply (function (x) x$correlations[c(3:30,32:34), c(3:30,32:34)] ) %>% ldply (function (x) ic.max.sim(x) )

dimnames(mean.compare.cov)[[1]] <- c(paste("Mean", mean.compare.cov[,1], sep = "."))
dimnames(min.compare.cov)[[1]] <- c(paste("Min", min.compare.cov[,1], sep = "."))
dimnames(max.compare.cov)[[1]] <- c(paste("Max", max.compare.cov[,1], sep = "."))
dimnames(ic.compare.cov)[[1]][c(1,3)] <- c(paste("ic.min", ic.compare.cov[,1][c(1,3)], sep = "."))
dimnames(ic.compare.cov)[[1]][c(2,4)] <- c(paste("ic.max", ic.compare.cov[,1][c(2,4)], sep = "."))
dimnames(ic.max.compare.cov)[[1]] <- c(paste("Max", ic.max.compare.cov[,1], sep = "."))
dimnames(ic.min.compare.cov)[[1]] <- c(paste("Min", ic.min.compare.cov[,1], sep = "."))

range.compare.large <- rbind(mean.compare.cov, ic.min.compare.cov, ic.max.compare.cov, ic.compare.cov)
range.compare.large <- range.compare.large[,-1]
range.compare.large <- as.data.frame(t(range.compare.large))


sample.size.list <- n.size[mask,2][c(3:30,32:34)]
range.compare.large$n <- sample.size.list
range.compare.large$.id <- rownames(range.compare.large)

str(range.compare.large)

#range.compare.backup <- range.compare
range.compare <- range.compare.large
range.compare <- range.compare.backup
#sp.main.data$Microcebus_griseorufus$rarefaction$rs
  plot.df <- melt(temp.rar.rs , value.name = "Mean.BS.RS")
  names(plot.df) <- c("Mean.BS.RS", "n") 
  plot.df <- as.data.frame(lapply(plot.df, as.numeric))
Plot.Raref.rs<- 
    ggplot(data = plot.df, aes(x= n, y = Mean.BS.RS)) + 
    geom_boxplot(aes( group = n), color="grey", outlier.shape = 8 ) + 
    geom_point(data = range.compare, aes(x = n, y = Mean.BS.RS), shape=17, size = 2.2, color = "red") +
    #geom_text(data = range.compare, aes(x = n, y = Mean.BS.RS, label = .id), color = "blue", size = 3, alpha = 0.3, angle = 70, hjust =  0.6) +
    geom_errorbar(data = range.compare, aes(x= n, ymin = Min.BS.RS, ymax = Max.BS.RS), alpha = 0.3, size = 1.2, color = "red") +
    geom_hline(data = range.compare, aes(yintercept = mean (Mean.BS.RS), y = mean (Mean.BS.RS)), alpha = 0.1, size = 2.5, color = "red") +
   ggtitle("Random Skewers") +
  theme_bw() +
    theme(plot.title = element_text(lineheight=.8, face="bold"),
            axis.text.x = element_text(size = 15), 
            axis.text.y = element_text(size = 15)) +
    ylab(label = "") + 
  xlab("Number of sampled specimens") +
   #xlab("Number of sampled specimens") + 
    scale_y_continuous(limits = c(0,1), breaks = c(0.00, 0.25, 0.5, 0.66, 0.75, 1.00)) +
    scale_x_continuous(limits = c(0, 68)) 
   
    

    plot.df.krz <- melt(temp.rar.krz, value.name = "Mean.BS.KRZ")
    names(plot.df.krz) <- c("Mean.BS.KRZ", "n") 
    plot.df.krz <- as.data.frame(lapply(plot.df.krz, as.numeric))
Plot.Raref.krz<- 
    ggplot(data = plot.df.krz, aes(x= n, y = Mean.BS.KRZ)) + 
      geom_boxplot(aes( group = n), color="grey", outlier.colour = "grey", outlier.shape = 8) + 
      geom_point(data = range.compare, aes(x = n, y = Mean.BS.KRZ), shape=17, size = 2.2, color = "blue") +
      #geom_text(data = range.compare, aes(x = n, y = Mean.BS.KRZ, label = .id), color = "blue", size = 3, alpha = 0.3, angle = 70, hjust =  0.6) +
      geom_errorbar(data = range.compare, aes(x= n, ymin = Min.BS.KRZ, ymax = Max.BS.KRZ), alpha = 0.3, size = 1.2, color = "blue") +
      geom_hline(data = range.compare, aes(yintercept = mean (Mean.BS.KRZ), y = mean (Mean.BS.KRZ)), alpha = 0.1, size = 2.5, color = "blue") +
   ggtitle("Krzanowski") +
  theme_bw() +    
  theme(plot.title = element_text(lineheight=.8, face="bold"),
            axis.text.x = element_text(size = 15), 
            axis.text.y = element_text(size = 15)) +
      ylab("") + 
  #  xlab("") +
  xlab("Number of sampled specimens") + 
      scale_y_continuous(limits = c(0,1), breaks = c(0.00, 0.25, 0.5, 0.75, 0.82, 1.00)) +
      scale_x_continuous(limits = c(0, 68)) 
      
Plot.Raref <- plot_grid(Plot.Raref.rs, 
                               #rarefaction.plot.rs.large, 
                               Plot.raref.krz, 
                               #rarefaction.plot.krz.large, 
                               labels = LETTERS[1:2], ncol=2)

#Antigo plot de rarefaçao
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
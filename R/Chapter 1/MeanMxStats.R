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
MMxStats$.sp %<>% gsub("_", ' ', .) 
MMxStats[,1] <- factor(unique(MMxStats[,1]), levels = unique(MMxStats[,1])[order(unique(MMxStats[,1]), decreasing = T)])
names(MMxStats)[2] <- "R^2"
names(MMxStats)[3] <- "PC1.percent"

MMxStats %>% #.[, c(1:3,9)] %>% 
  gather(key= .MMxStats, value=value, 2:3) %>% 
  ggplot(aes(y=  value, x = .sp, fill = .MMxStats ), varwidth = F) +
  geom_bar(aes( y=  value, x = .sp, fill = .MMxStats ), stat = "identity", position="dodge")+
  theme_bw() +
  #ggtitle("Mean matrix evolutionary statistics by specie") +
  theme(axis.text.y = element_text(face = "italic")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  scale_fill_grey() +
  #scale_fill_brewer(palette="Spectral") +
  scale_colour_brewer(name = "Mean Matrix Statistics") + 
  facet_wrap(~.MMxStats) +
  coord_flip() +
  theme(legend.position="none")


MMxStats<- cov.mx[mask] %>% ldply(function(x) MeanMatrixStatistics(x))
names(MMxStats)[2] <- "MeanSquaredCorrelation"
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

gm.mean <- sp.main.data[mask] %>% ldply(function(x) x$gm.mean)

MMxStats.gm<- cbind(MMxStats, gm.mean[mask,2])
names(MMxStats.gm)[11] <- "Geometric.mean"

names (MMxStats.gm)[3] <- "PC1.percent"
#################################################################################################################
####################### ESSE AQUI VAI ENTRAR: NENHUMA RELAÇAO ENTRE TAMANHO DO CRANIO E PC1% ####################
#################################################################################################################
Plot.GM.PC1 <- MMxStats.gm %>%
  ggplot( ., aes(x = Geometric.mean, y = PC1.percent, label = .sp), varwidth = T) +
  geom_text(size =2, vjust = 1, alpha = 0.8)  +
  stat_smooth(method="lm", aes(group=1), color = "red") +
  scale_y_continuous(limits = c(0.15, 0.77), breaks = c(0.25, 0.5, 0.75) )+
  theme_bw() +
  geom_point(pch = 20) +
  geom_text(x = 20,  y = 0.64, label = "r.squared = 0.001853 \n 40 DF,  p-value: 0.7867", color = "red", size = 2) + 
  #ggtitle("Geometric mean x PC1.percent") +
  theme(legend.position  ="none") +
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  xlab("Geometric mean") +
  ylab("PC1 % of variance")

summary(lm(Geometric.mean ~ PC1.percent , data = MMxStats.gm))

Plot.R2.PC1 <- MMxStats.gm %>%
  ggplot( ., aes(x = MeanSquaredCorrelation, y = PC1.percent, label = .sp), varwidth = T) +
  geom_text(size =2, vjust = 1, alpha = 0.8)  +
  stat_smooth(method="lm", aes(group=1), color = "red") +
  scale_y_continuous(limits = c(0.15, 0.77), breaks = c(0.25, 0.5, 0.75) )+
  theme_bw() +
  geom_point(pch = 20) +
  geom_text(x = 0.11,  y = 0.74, label = "r.squared = 0.8464 \n 40 DF,  p-value: < 2e-16", color = "red", size = 2) + 
  #ggtitle("Geometric mean x PC1.percent") +
  theme(legend.position  ="none") +
  theme(plot.title = element_text(lineheight=.8, face="bold"))+
  xlab("Mean squared correlation (r²)") +
  ylab("PC1 % of variance")

summary(lm(MeanSquaredCorrelation ~ PC1.percent , data = MMxStats.gm))

EigVal <- cov.mx[mask] %>% ldply(function(x) eigen (x)$ values) 
MMxStats.gm$EigVal1 <- EigVal$V1

Plot.GM.EigenVal1 <- MMxStats.gm %>%
  ggplot( ., aes(x = Geometric.mean, y = EigVal1, label = .sp), varwidth = T) +
  geom_text(size =2, vjust = 1, alpha = 0.8)  +
  stat_smooth(method="lm", aes(group=1), color = "red") +
  #scale_y_continuous(limits = c(0.15, 0.77), breaks = c(0.25, 0.5, 0.75) )+
  theme_bw() +
  geom_point(pch = 20) +
  geom_text(x = 20,  y = 0.64, label = "r.squared = 0.1986 \n 40 DF,  p-value: 0.0031", color = "red", size = 2) + 
  #ggtitle("Geometric mean x PC1.percent") +
  theme(legend.position  ="none") +
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  xlab("Geometric mean") +
  ylab("PC1")

summary(lm(Geometric.mean ~ EigVal1 , data = MMxStats.gm))
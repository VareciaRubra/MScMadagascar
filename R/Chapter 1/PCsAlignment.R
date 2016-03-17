########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
current.data <- sp.main.data
PCs1to4<- current.data[mask] %>% llply(function(x) as.data.frame(eigen(x$matrix$cov)$vectors[,1:4]) ) %>% llply(function(x) as.list(x) )
PC1 <- current.data[mask] %>% ldply(function(x) eigen(x$matrix$cov)$vectors[,1])
names.sp <-  PC1[,1]
PC1<- as.matrix(PC1[,-1])
dimnames(PC1)[[1]] <- names.sp
# se estiver correlacionado negativamente com o vetor isometrico, multiplicar por -1
ed.names <- names(current.data[mask][[1]]$ed)
Iso.Compare <- function(x) {
  isometrico<- rep( (1/sqrt(39)) , 39) # construindo vetor isométrico 
  corr = rep(NA, length(x) )
  re.oriented = x
  for (i in 1:length(x)) {
    corr[i]<- x[[i]] %*% isometrico
    
    if(corr[i] <0) {re.oriented[[i]] <- as.numeric(as.matrix(x[[i]]) * -1)}  else x[[i]]
  }
  
  names (corr) <- paste ("PC", 1:length(x), sep = "")
  names (re.oriented) <- paste ("PC", 1:length(x), sep = "") 
  
  return(results= list("corr" = corr, "re.oriented" = re.oriented)) 
}
PCs1to4<- current.data[mask] %>% llply(function(x) as.data.frame(eigen(x$matrix$cov)$vectors[,1:4]) ) %>% llply(function(x) as.list(x) )
PCs1to4.log<- current.data[mask] %>% llply(function(x) as.data.frame(eigen(x$matrix$cov.log)$vectors[,1:4]) ) %>% llply(function(x) as.list(x) ) 
### com os valores de PCs Logados
myPalette <- colorRampPalette(rev(brewer.pal(5, 'Spectral')), space = 'Lab')(n = 1000)

Iso.Compare.reoriented <- llply(PCs1to4, Iso.Compare) %>% ldply(function(x) as.data.frame(x$re.oriented)) 
names(Iso.Compare.reoriented) <- c(".sp", "PC1", "PC2", "PC3", "PC4")
Iso.Compare.reoriented$.ed <- ed.names
Iso.Compare.reoriented$.sp <- factor(Iso.Compare.reoriented$.sp, levels = unique(Iso.Compare.reoriented$.sp)) 
Iso.Compare.reoriented$.ed <- factor(Iso.Compare.reoriented$.ed, levels = rev(ed.names) ) 
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")

Iso.Compare.cor.iso <- llply(PCs1to4, Iso.Compare) %>% ldply(function(x) x$corr) 
names(Iso.Compare.cor.iso) <- c(".sp", "PC1", "PC2", "PC3", "PC4")
Iso.Compare.cor.iso$.sp %<>% gsub("_", " ", .)
Iso.Compare.cor.iso$.sp <- factor(Iso.Compare.cor.iso$.sp, levels = rev(unique(Iso.Compare.cor.iso$.sp))) 
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
########################################################################################################
########################################################################################################
########################################################################################################
Iso.Compare.cor.iso %>% gather(key="Isometric.Correlation", value=value, 2:5 ) %>%
  ggplot (.) +
  geom_tile(aes(x = Isometric.Correlation, y = .sp, fill = abs(value) ), alpha = 0.9 ) +
  #geom_text(aes(x= Isometric.Correlation, y = .sp, label = abs(round(value, digits = 2) ) ), color = "black" ) +
  theme_bw() +
  scale_fill_gradientn(name = 'Correlation with \nIsometric vector', colours = myPalette(50)) +
  ylab ('') + xlab ('') +
  #scale_x_discrete(limits = levels(PCs1to4$.pcScore), breaks = c(0.2, 0.4, 0.6, 0.8, 1) ) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 8),
        axis.ticks = element_line(size = 0),
        axis.text.y = element_text(face = "italic", size= 10),
        plot.title = element_text(lineheight=.8, face="bold")) 

Iso.Compare.reoriented %>% gather(key="Isometric.Correlation", value=value, 2:5 ) %>%
  ggplot (.) +
  geom_tile(aes(x = Isometric.Correlation, y = .sp, fill = abs(value) ) ) +
  geom_text(aes(x= Isometric.Correlation, y = .sp, label = round(value, digits = 2 ) ) ) +
  theme_bw() +
  scale_fill_gradientn(name = 'Correlation with \nIsometric vector', colours = myPalette(50)) +
  ylab ('') + xlab ('') +
  #scale_x_discrete(limits = levels(PCs1to4$.pcScore), breaks = c(0.2, 0.4, 0.6, 0.8, 1) ) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10),
        axis.ticks = element_line(size = 0),
        axis.text.y = element_text(face = "italic", size= 10),
        plot.title = element_text(lineheight=.8, face="bold")) 



size.matters <- sp.main.data[mask] %>% ldply(., function(x) x$rs.size.comparisson[[1]][lower.tri(x$rs.size.comparisson[[1]])] ) 
colnames(size.matters) <- c(".sp", "rawXsizeless", "rawXlog", "sizelessXlog")
size.matters <- cbind(size.matters, Iso.Compare.cor.iso)


size.matters <- cbind(size.matters %>% melt(id.vars = 1, measure.vars = c(2:4), variable_name = "size.compare") , abs(rep(Iso.Compare.cor.iso[,2], 3) ) )
colnames(size.matters)[4] <- "Pc.iso.cor"
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")

size.matters %>% 
  ggplot(  aes(x= variable, y = value) ) + 
  geom_violin() + 
  geom_text(aes(label = .sp, color = Pc.iso.cor )  ) +
  scale_color_gradientn(name = 'Pc1 Correlation with \nIsometric vector', colours = myPalette(50)) 

size.matters %>% 
  ggplot(  aes(y= variable, x = value) ) + 
  #geom_violin() + 
  geom_text(aes(y= Pc.iso.cor, x = value, label = .sp )  ) +
  #scale_color_gradientn(name = 'Pc1 Correlation with \nIsometric vector', colours = myPalette(50)) 
  facet_wrap(~variable)
  

scale_color_gradient(name = 'PC1 Correlation with \nIsometric vector',low = myPalette[1], high = myPalette[50], mid = myPalette[25]) 

size.matters %>%
  geom_boxplot(x )



sp.main.data[mask] %>% llply(., function(x) x$rs.size.comparisson[[1]] ) 





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

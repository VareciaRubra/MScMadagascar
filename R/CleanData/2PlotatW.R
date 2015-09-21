#Entrando a árvore
#Arvore
treefile = read.nexus(file = "attaches/fbd421agerange_edited.tre")
species <- treefile$tip.label[treefile$tip.label %in% names(sample.no.na)]
pruned.tree<-drop.tip(treefile,treefile$tip.label[-match(species, treefile$tip.label)])
plot(pruned.tree)
nodelabels()
#matrizes ancestrais ponderadas pela filogenia para cada nó da arvore: 
#fazer a chamada pelo numero correspondente do nó na filogenia
Ancestral.Matrices<- PhyloW(tree = pruned.tree, tip.data = cov.no.na, tip.sample.size = sample.no.na)
#matriz W selecionada a partir da PhyloW do pacote ATENÇÃO QUE A GRAFIA TEM QUE SER COM ASPAS PARA NUMERO DO NÓ
Wmat <- Ancestral.Matrices$"42" 
GeralMorphoSpace = CalculateMatrix(manova(as.matrix(all.main.data$All$ed)  ~ Especie, data = as.data.frame(all.main.data$All$info) ) )
Wmat <- GeralMorphoSpace
#definindo qual o conjunto de dados que quero plotar no espaçõ da matriz determinada:
current.data <- all.main.data #todo o banco de dados
# current.data <- extant.main.data #apenas os viventes
# current.data<- madagascar.main.data # apenas os de madagascar
# current.data <- extant.madagascar.main.data # apenas os viventes de madagascar

current.data$All$info$Familia <- factor(current.data$All$info$Familia, levels = unique(current.data$All$info$Familia))

#juntando as informações de cada indivíduo com suas respectivas distancias
lin_data = ldply(current.data, function(x) tbl_df(cbind(x$info, x$ed))) 
#Complete info
lin_data %<>% mutate (., .Cinfo = paste(Status, Regiao, Familia, Genero, Especie, sep = ".")) #criando uma coluna de info que vai ter dados taxonomicos de familia e genero (para depois separar o plot em cores)
#Resumed info
lin_data %<>% mutate (., .Rinfo = paste(Status, Regiao, Familia, Genero, sep = ".")) #criando uma coluna de info que vai ter dados taxonomicos de familia e genero (para depois separar o plot em cores)
#myformula = paste0('~', paste(names(select(lin_data, IS_PM:BA_OPI)), collapse = '+'))
rownames(lin_data) <- lin_data$Tombo

######## Biplot dos scores de cada indivíduo no espaço da W desejada##############
lin_data <- select(lin_data, .Cinfo, .Rinfo, Tombo, Especie, IS_PM:BA_OPI)
#resp <- as.data.frame(as.matrix(na.omit(select(current.data$All$ed, IS_PM:BA_OPI) ) ) %*% eigen(Wmat)$vectors[,c(1:39)])
mx.all <- as.matrix(na.omit(select(current.data$All$ed, IS_PM:BA_OPI) ) )
mx.all<- scale(mx.all, center = TRUE, scale = FALSE)
resp <- as.data.frame(mx.all %*% eigen(Wmat)$vectors[,c(1:39)])
resp %<>% mutate(., .id = rownames(resp) )
names(resp) %<>% gsub("V", "PC", .)
resp.info <- select(lin_data, c(.Cinfo, .Rinfo, Tombo) )
resp.info <- resp.info[resp.info$Tombo %in% resp$.id, ]
plot.W <- cbind(resp.info, resp)
plot.W$PC1 <- plot.W$PC1 *(-1) + 120
plot.W$PC2 <- plot.W$PC2 *(-1)
plot.W$PC4 <- plot.W$PC4 *(-1)

plot.W.spread <- ddply(plot.W, .(.Rinfo), plyr::summarise, "hPC1"=PC1[chull(PC1,PC4)],
      "hPC4"=PC4[chull(PC1,PC4)])
plot.W.spread %<>% separate(.Rinfo, c('Status', 'Regiao', 'Family', 'Genero'), sep = "\\.")
hulls <- plot.W.spread
hulls$Family <- factor(hulls$Family, levels = unique(current.data$All$info$Familia))

hulls <- ddply(plot.W, .(.Cinfo), plyr::summarise, "hPC1"=PC1[chull(PC1,PC4)],
               "hPC4"=PC4[chull(PC1,PC4)])
hulls %<>% separate(.Cinfo, c('Status', 'Regiao', 'Family', 'Genero', 'Especie'), sep = "\\.")
hulls$Family <- factor(hulls$Family, levels = unique(current.data$All$info$Familia))

points.info <- ddply(plot.W, .(.Cinfo), numcolwise(mean))
points.info %<>% separate(.Cinfo, c('Status', 'Regiao', 'Family', 'Genero', 'Especie'), sep = "\\.")


pc_plot <- ggplot(plot.W, aes(PC1, PC4)) +
  geom_polygon(aes(hPC1, hPC4, fill = Family, color = Family, group= Especie, log = "x"), 
               data = hulls, alpha=0.2) + 
  geom_point(data = points.info,
             aes(PC1, PC4, group= Especie, shape = interaction(Status, Regiao), color = Family), 
             size = 5, alpha = 0.3, log = "x") + 
  scale_shape(name = "Localization and status", labels = c("Madagascar - Extants", "Madagascar - Extincts", "Out Madagascar -  Extants")) +
  geom_text(data = ddply(points.info, .(Genero), numcolwise(mean)),
            aes(PC1, PC4, label= Genero, size = log((abs(PC1) + 200)) ),  alpha = 0.4) +
#   scale_size(name = "Size reference\nAverage PC1 score per genus", breaks = c(5.2, 5.6, 6.0),
#              labels = c("Small", "Medium", "Large") ) +
  theme_bw() + 
#   theme(axis.title.x = element_text(face="bold", size=10) ) + labs(x = "PC1 = Size") +
#   theme(axis.title.y = element_text(face="bold", size=10) ) + labs(y = "PC2 = Snout Lenght x Vault Volume") +
  ggtitle("Projection at Ancestral Matrix's Morphospace") +
  #coord_cartesian(ylim=c(-63, 63), xlim = c(-10, 550)) +
  theme(plot.title = element_text(lineheight=.8, face="bold")) 
  pc_plot + coord_trans(x = "log10") 
  
  
  
  
  
  plot.W.spread <- ddply(plot.W, .(.Rinfo), plyr::summarise, "hPC1"=PC1[chull(PC1,PC2)],
                         "hPC2"=PC2[chull(PC1,PC2)])
  plot.W.spread %<>% separate(.Rinfo, c('Status', 'Regiao', 'Family', 'Genero'), sep = "\\.")
  hulls <- plot.W.spread
  hulls$Family <- factor(hulls$Family, levels = unique(current.data$All$info$Familia))
  
  hulls <- ddply(plot.W, .(.Cinfo), plyr::summarise, "hPC1"=PC1[chull(PC1,PC2)],
                 "hPC2"=PC2[chull(PC1,PC2)])
  hulls %<>% separate(.Cinfo, c('Status', 'Regiao', 'Family', 'Genero', 'Especie'), sep = "\\.")
  hulls$Family <- factor(hulls$Family, levels = unique(current.data$All$info$Familia))
  
  points.info <- ddply(plot.W, .(.Cinfo), numcolwise(mean))
  points.info %<>% separate(.Cinfo, c('Status', 'Regiao', 'Family', 'Genero', 'Especie'), sep = "\\.")
  
  
  pc_plot <- ggplot(plot.W, aes(PC1, PC2)) +
    geom_polygon(aes(hPC1, hPC2, fill = Familia, color = Familia, group= Especie, log = "x"), 
                 data = hulls, alpha=0.2) + 
    geom_point(data = points.info,
               aes(PC1, PC2, group= Especie, shape = interaction(Status, Regiao), color = Familia), 
               size = 5, alpha = 0.3, log = "x") + 
    scale_shape(name = "Localization and status", labels = c("Madagascar - Extants", "Madagascar - Extincts", "Out Madagascar -  Extants")) +
    geom_text(data = ddply(points.info, .(Genero), numcolwise(mean)),
              aes(PC1, PC2, label= Genero, size = log((abs(PC1) + 200)) ),  alpha = 0.4) +
    #   scale_size(name = "Size reference\nAverage PC1 score per genus", breaks = c(5.2, 5.6, 6.0),
    #              labels = c("Small", "Medium", "Large") ) +
    theme_bw() + 
    #   theme(axis.title.x = element_text(face="bold", size=10) ) + labs(x = "PC1 = Size") +
    #   theme(axis.title.y = element_text(face="bold", size=10) ) + labs(y = "PC2 = Snout Lenght x Vault Volume") +
    ggtitle("Projection at Ancestral Matrix's Morphospace") +
    #coord_cartesian(ylim=c(-63, 63), xlim = c(-10, 550)) +
    theme(plot.title = element_text(lineheight=.8, face="bold")) 
  pc_plot + coord_trans(x = "log10") 
  
  


#################################################################################################################
# INTERPRETANDO OS PCs DA W
#################################################################################################################
#Escolha a matriz que voce quer, meu bem:
plot(pruned.tree)
nodelabels()
Wmat <- Ancestral.Matrices$"42" 
Wmat <- GeralMorphoSpace
WPCs <- as.matrix(eigen(Wmat)$vectors)
dimnames(WPCs)[[2]] <- paste("PC", 1:39, sep= "")
dimnames(WPCs)[[1]] <- dimnames(all.main.data$All$matrix$cov)[[1]]

m.WPCs <- melt(WPCs)
m.WPCs$Var1<- factor( m.WPCs$Var1, levels = levels( m.WPCs$Var1)[39:1])
m.WPCs.position =  m.WPCs
m.WPCs.position$Var1 <- as.numeric( m.WPCs.position$Var1)
m.WPCs.position$Var2 <- as.numeric( m.WPCs.position$Var2)
m.WPCs.position$value= round( m.WPCs.position$value, 3)

myPalette <- colorRampPalette(brewer.pal(11, 'Spectral')[-(4:6)], space = 'Lab')(n = 10)
ggplot ( m.WPCs) +
  geom_tile(aes(x = Var2, y = Var1, fill = value)) +
  scale_fill_gradientn(name = '', colours = myPalette) +
  ylab ('') + xlab ('') + labs(title = "GeralMorphoSpace PCs com fosseis") + 
  geom_text(data =  m.WPCs.position, size = 3, aes(x = Var2, y = Var1, label = round(value, 2)) ) + 
  theme(axis.text.x = element_text(angle = 90),
        #axis.text.y = element_blank(),
        axis.ticks = element_line(size = 0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        rect = element_blank(), line = element_blank())
#################################################################################################################

KrzCor(cov.x = as.matrix(eigen(Ancestral.Matrices$"42")$vectors) , cov.y = as.matrix(eigen(GeralMorphoSpace)$vectors) )

#Entrando a árvore
#Arvore
treefile = read.nexus(file = "attaches/fbd421agerange_edited.tre")
species <- treefile$tip.label[treefile$tip.label %in% names(sample.no.na)]
pruned.tree<-drop.tip(treefile,treefile$tip.label[-match(species, treefile$tip.label)])
plot(pruned.tree)
nodelabels()
#matrizes ancestrais ponderadas pela filogenia para cada nó da arvore: 
#fazer a chamada pelo numero correspondente do nó na filogenia
#Ancestral.Matrices<- PhyloW(tree = pruned.tree, tip.data = cov.no.na, tip.sample.size = sample.no.na)
#matriz W selecionada a partir da PhyloW do pacote ATENÇÃO QUE A GRAFIA TEM QUE SER COM ASPAS PARA NUMERO DO NÓ
Wmat <- Gen.cov.list$W.Prosimian
#all.main.data$All$info [aaply (as.matrix (all.main.data$All$ed), 1, function (L) any (is.na(L))), ]
GeralMorphoSpace = CalculateMatrix(manova(as.matrix(all.main.data$All$ed) ~ Especie, data = as.data.frame(all.main.data$All$info) ))
Wmat <- GeralMorphoSpace
#definindo qual o conjunto de dados que quero plotar no espaçõ da matriz determinada:
current.data <- all.main.data #todo o banco de dados

current.data$All$info$Familia <- factor(current.data$All$info$Familia, levels = unique(current.data$All$info$Familia))

#juntando as informações de cada indivíduo com suas respectivas distancias
lin_data = ldply(current.data, function(x) tbl_df(cbind(x$info, x$ed))) 
lin_data <- lin_data %>% filter(., Status != "Extinct")

#Complete info
lin_data %<>% mutate (., .Cinfo = paste(Status, Regiao, Familia, Genero, Especie, sep = ".")) #criando uma coluna de info que vai ter dados taxonomicos de familia e genero (para depois separar o plot em cores)
#Resumed info
lin_data %<>% mutate (., .Rinfo = paste(Status, Regiao, Familia, Genero, sep = ".")) #criando uma coluna de info que vai ter dados taxonomicos de familia e genero (para depois separar o plot em cores)
#myformula = paste0('~', paste(names(select(lin_data, IS_PM:BA_OPI)), collapse = '+'))
lin_data$Tombo <- factor(x = lin_data$Tombo, levels = unique(lin_data$Tombo))
rownames(lin_data) <- lin_data$Tombo

######## Biplot dos scores de cada indivíduo no espaço da W desejada##############
lin_data <- select(lin_data, .Cinfo, .Rinfo, Tombo, Especie, IS_PM:BA_OPI)
mx.all <- as.matrix(na.omit(select(current.data$All$ed[current.data$All$info$Status !=  "Extinct",], IS_PM:BA_OPI) ) )
mx.all<- scale(mx.all, center = TRUE, scale = FALSE)
resp <- as.data.frame(mx.all %*% eigen(Wmat)$vectors[,c(1:39)])
resp %<>% mutate(., .id = rownames(resp) )
names(resp) %<>% gsub("V", "PC", .)
resp.info <- select(lin_data, c(.Cinfo, .Rinfo, Tombo) )
resp.info <- resp.info[resp.info$Tombo %in% resp$.id, ]
table(as.character(resp.info$Tombo) == resp$.id)
plot.W <- cbind(resp.info, resp)
plot.W$PC1 <- plot.W$PC1 *(-1) + 90
plot.W$PC2 <- plot.W$PC2 *(-1) 
plot.W$PC3 <- plot.W$PC3 *(-1) 
plot.W$PC4 <- plot.W$PC4 *(-1) + 30

  hulls <- ddply(plot.W, .(.Cinfo), plyr::summarise, "hPC1"=PC1[chull(PC1,PC4)],
                 "hPC4"=PC4[chull(PC1,PC4)], "tombo"=Tombo[chull(PC1,PC4)])
  
  hulls <- ddply(plot.W, .(.Cinfo), plyr::summarise, "hPC1"=PC1[chull(PC1,PC4)],
                 "hPC4"=PC4[chull(PC1,PC4)])
  hulls %<>% separate(.Cinfo, c('Status', 'Regiao', 'Family', 'Genero', 'Especie'), sep = "\\.")
  hulls$Family <- factor(hulls$Family, levels = unique(current.data$All$info$Familia))
  
  points.info <- ddply(plot.W, .(.Cinfo), numcolwise(mean))
  points.info %<>% separate(.Cinfo, c('Status', 'Regiao', 'Family', 'Genero', 'Especie'), sep = "\\.")
  
  pc_plot14 <- ggplot(plot.W, aes(PC1, PC4)) +
    #geom_polygon(aes(hPC1, hPC4, fill = Family, color = Family, group= Especie, log = "x"), data = hulls, alpha=0.2) + 
    geom_point(data = points.info,
               aes (PC1, PC4, group= Genero, shape = Regiao, color = Family, log = "x"), 
               alpha = 0.7, size = 2) + 
    #theme(legend.text = element_text(angle = 90, size = 10, hjust = 3, vjust = 3, face = 'italics')) +
    #scale_shape(name = "Localização e status", labels = c("Madagascar - viventes", "Madagascar - Extintos", "Fora de Madagascar - Viventes")) +
    geom_text(data = ddply(points.info, .(Genero, Family), numcolwise(mean)),
              aes(PC1, PC4, label= Genero, color = Family),  alpha = 0.4 ) +
#        scale_size(name = "Size reference\nAverage PC1 score per genus", breaks = c(5.2, 5.6, 6.0),
#                labels = c("Small", "Medium", "Large") ) +
    scale_colour_brewer( palette="Set1") +
    theme_bw() + 
    #labs(x = "PC1 = Tamanho") +
    #labs(y = "PC4 = Contraste entre comprimento\nda face e volume do neurocrânio ") +
    ggtitle("Projeção no Morfoespaço da Matriz Ancestral") +
    #coord_cartesian(ylim=c(-2, 58), xlim = c(-10, 580)) +
    coord_fixed() +
    guides(size = FALSE, 
           text = FALSE, 
           colour = guide_legend(title.position = "top", title = "Família", size = 40, override.aes = list(alpha = 1) ), 
           fill = guide_legend(title.position = "top", title = "Família", size = 40, override.aes = list(alpha = 1) ),
           shape = guide_legend(title.position = "top", override.aes = list(alpha = 1, size = 5), nrow = 1) ) +
    theme(legend.position = "bottom", 
          axis.title= element_text(size=15),
          plot.title = element_text(lineheight=.8, size = 13))

  plot_grid(plot_grid(pc_plot12+ theme(legend.position = "none"), 
            pc_plot13+ theme(legend.position = "none"), 
            pc_plot23+ theme(legend.position = "none"), ncol =1),
            pc_plot14, ncol = 2)
  
  pc_plot  + coord_trans(x = "log10", y = "log10")
  pc_plot + coord_fixed()



  
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

KrzCor(cov.x = Ancestral.Matrices$"42" , cov.y = GeralMorphoSpace )
RandomSkewers(cov.x = Ancestral.Matrices$"42" , cov.y = GeralMorphoSpace )



## Canonical variates
ProjetaDados = function(y,var.y){
  y = as.matrix(y)
  n = nrow(y)
  p = ncol(y)
  eigen.y = eigen(var.y)
  eVal = eigen.y$values
  eVec = eigen.y$vectors
  Scores = array(0., c(n,p))
  for( i in 1:n){
    Scores[i,] = t(eVec)%*%(as.numeric(y[i,]))
  }
  for(i in 1:p){
    Scores[,i] <- Scores[,i]/eVal[i]
  }
  return(Scores)
}

data_selected <- lin_data %>% 
  select(., .Cinfo, Tombo, IS_PM:BA_OPI) %>% 
  separate(., .Cinfo, c('Status', 'Regiao', 'Family', 'Genero', 'Especie'), sep = "\\." ) %>% 
  filter ( Family == "Lemuridae") %>%
  filter ( Genero == "Eulemur") %>%
  mutate (., .Cinfo = paste(Status, Regiao, Family, Genero, Especie, sep = ".")) 
data_selected$Tombo <- factor(x =data_selected$Tombo, levels = unique(lin_data$Tombo))
rownames(data_selected) <- data_selected$Tombo

#Escolha a matriz que voce quer, meu bem:
plot(pruned.tree)
nodelabels()
Wmat <- Ancestral.Matrices$"65" 

lm.within = lm(as.matrix(select(data_selected, IS_PM:BA_OPI)) ~ select(data_selected, Especie)[,1])

current.data_projected_W = cbind(select(data_selected, .Cinfo, Tombo),  
                                 ProjetaDados(select(data_selected, IS_PM:BA_OPI), Wmat))

Bmat = cov(daply(na.omit(current.data_projected_W), .(.Cinfo), function(x) colMeans(x[,-c(1, 2)])))

resp <- cbind(select(data_selected, .Cinfo, Especie, Tombo),  scale(ProjetaDados(select(data_selected, IS_PM:BA_OPI), Wmat) %*% eigen(Bmat)$vectors[,1:4], scale = TRUE))
resp <- resp[!is.na(resp[,3]),]
names(resp) <- c(".Cinfo", ".sp", "Tombo", "CV1", "CV2", "CV3", "CV4")
table(rownames(resp) == resp$Tombo) 
resp %<>% mutate(., .Cinfo = paste (.Cinfo, Tombo, sep = ".") )

hulls <-ddply(resp, .(.Cinfo), plyr::summarise, "hpc1"=CV1[chull(CV1,CV2)],
              "hpc2"=CV2[chull(CV1,CV2)])
hulls %<>% separate(.Cinfo, c('Status', 'Regiao', 'Family', 'Genero', 'Especie', 'Tombo'), sep = "\\.")
rownames(hulls) <- hulls$Tombo
resp %<>% ddply(., .(.Cinfo), numcolwise(mean))
resp %<>% separate(.Cinfo, c('Status', 'Regiao', 'Family', 'Genero', 'Especie', 'Tombo'), sep = "\\.")

cv_plot_12 <- ggplot(resp, aes(CV1, CV2)) +
  geom_polygon(aes(hpc1, hpc2, fill = Especie, group= Especie), hulls, alpha=.3) + 
  geom_point(data = resp,
             aes(CV1, CV2, group= Especie, color = Especie, shape = Genero, size= Status ), size = 5, alpha = 0.5) +
  geom_text(data = hulls,
             aes(hpc1, hpc2, group= Especie, label = Tombo), size = 5, alpha = 0.5) +
  theme_bw() +
  ggtitle("Canonical Variates") +
  theme(plot.title = element_text(lineheight=.8, face="bold")) 
cv_plot_12

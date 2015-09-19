####### Selecionando o que que vc quer plotar #########
current.data <- madagascar.main.data

lin_data = ldply(current.data, function(x) tbl_df(cbind(x$info, x$ed))) 
lin_data %<>% mutate (., .info = paste(Familia, Genero, Especie, sep = "."))
individual_ID <- unlist(llply(current.data, function(x) x$info$Tombo))
myformula = paste0('~', paste(names(select(lin_data, IS_PM:BA_OPI)), collapse = '+'))
rownames(lin_data) <- individual_ID
PRCOMP <- princomp(as.formula(myformula), data = lin_data, scale = T)
#PRCOMP$scores
######## Biplot PC1 x PC2 ##############
current.data <- select(lin_data, .info, IS_PM:BA_OPI)
rownames(current.data) <- rownames(lin_data)
#PRCOMP %>% biplot
#resp %<>% mutate(., PC1 = PRCOMP$scores[,1], PC2=PRCOMP$scores[,2]) 
#The diffent grids for ploting in W space: all, only extants, and only madagascar extant and with extinct.
Wmat.All.Fuckers = CalculateMatrix(manova(as.matrix(all.main.data$All$ed)  ~ Especie, data = as.data.frame(all.main.data$All$info) ) )
Wmat.extant.madagascar.Fuckers = CalculateMatrix(manova(as.matrix(extant.madagascar.main.data$All$ed)  ~ Especie, data = as.data.frame(extant.madagascar.main.data$All$info) ) )
Wmat = CalculateMatrix(manova(as.matrix(current.data$All$ed)  ~ Especie, data = as.data.frame(current.data$All$info) ) )
GeralMorphoSpace = cov(as.matrix(current.data$All$ed))


resp <- as.data.frame(as.matrix(na.omit(select(current.data, IS_PM:BA_OPI) ) ) %*% eigen(Wmat.extant.madagascar.Fuckers)$vectors[,c(1:39)])
resp %<>% mutate(., .id = rownames(resp) )
names(resp) %<>% gsub("V", "PC", .)
resp.info <- select(lin_data, c(.info,Arquivo:Data_dado) )
resp.info <- resp.info[resp.info$Tombo %in% resp$.id, ]
plot.W <- cbind(resp.info, resp)
hulls <- ddply(plot.W, .(.info), plyr::summarise, "hpc1"=PC1[chull(PC1,PC2)],
               "hpc2"=PC2[chull(PC1,PC2)])
hulls %<>% separate(.info, c('Familia', 'Genero', "Especie"), sep = "\\.")

# levels(plot.W$Especie) <- as.factor(unique(resp.info$Especie))
# levels(plot.W$Genero) <- as.factor(unique(resp.info$Genero))
# levels(plot.W$Familia) <- as.factor(unique(resp.info$Familia))
# hulls$Genero <- factor (hulls$Genero, levels = hulls$Genero )
hulls$Especie <- unique (hulls$Especie, levels = unique(as.factor(hulls$Especie) ) )
pc_plot <- ggplot(plot.W, aes(PC1, PC2)) +
  geom_polygon(aes(hpc1, hpc2, fill = Familia, color = Familia, group= Genero), data = hulls, alpha=.3) + 
  geom_point(data = ddply(plot.W, .(Especie), numcolwise(mean)),
             aes(PC1, PC2, group= Especie), size = 3) + 
  theme_bw() + ggtitle("Cranial traits Within-group PC scores") + 
  coord_fixed()
pc_plot


pc_plot <- ggplot(plot.W, aes(PC1, PC4))+
  geom_polygon(aes(hpc1, hpc2, fill = Familia, color = Familia, group= Especie ), data = hulls, alpha=.2) + 
  geom_point(data = ddply(plot.W, .(Tombo), numcolwise(mean)),
             aes(PC1, PC4, group= Tombo), size = 1, color = "grey", aplha = 0.3) +
  #   geom_text(data = ddply(plot.W, .(Tombo), numcolwise(mean)),
  #              aes(PC1, PC2, label= Tombo), size = 3) 
  #   geom_point(data = ddply(plot.W, .(Genero), numcolwise(mean)),
  #             aes(PC1, PC2, group= Genero), size = 3) +
  geom_text(data = ddply(plot.W, .(Genero), numcolwise(mean)),
            aes(PC1, PC4, label= Genero), size = 5) +
  theme_bw() + ggtitle("Cranial traits Within-group PC scores")
pc_plot


pc_plot <- ggplot(plot.W, aes(PC1, PC2)) +
  geom_polygon(aes(hpc1, hpc2, fill = Familia, color = Familia, group= Genero ), data = hulls, alpha=.2) + 
  geom_point(data = ddply(plot.W, .(Tombo), numcolwise(mean)),
             aes(PC1, PC2, group= Tombo), size = 1, color = "grey", aplha = 0.3) +
  #   geom_text(data = ddply(plot.W, .(Tombo), numcolwise(mean)),
  #              aes(PC1, PC2, label= Tombo), size = 3) 
  #   geom_point(data = ddply(plot.W, .(Genero), numcolwise(mean)),
  #             aes(PC1, PC2, group= Genero), size = 3) +
  geom_text(data = ddply(plot.W, .(Genero), numcolwise(mean)),
            aes(PC1, PC2, label= Genero), size = 5) +
  theme_bw() + ggtitle("Cranial traits Within-group PC scores")
pc_plot



Wmat.extant.Fuckers = CalculateMatrix(manova(as.matrix(extant.main.data$All$ed)  ~ Especie, data = as.data.frame(extant.main.data$All$info) ) )
Wmat.madagascar.Fuckers = CalculateMatrix(manova(as.matrix(madagascar.main.data$All$ed)  ~ Especie, data = as.data.frame(madagascar.main.data$All$info) ) )


# creating the masks to index the current set of matrices, means and repetabilities
mask.extant<- unique(All.raw.main.data$All$Especie[All.raw.main.data$All$Status != "Extinct"])
mask.madagascar<- unique(All.raw.main.data$All$Especie[All.raw.main.data$All$Regiao == "Madagascar"])
mask.madagascar.extant <- unique(All.raw.main.data$All$Especie[All.raw.main.data$All$Status != "Extinct" &All.raw.main.data$All$Regiao == "Madagascar"])

#selecionando todas as médias, dependendo do grid
means.All.Fuckers<- current.data %>% llply(function(x) x$ed.means )
means.extant.Fuckers<- means.All.Fuckers[mask.extant] 
means.madagascar.Fuckers<- means.All.Fuckers[mask.madagascar]
means.extant.madagascar.Fuckers<- means.All.Fuckers[mask.madagascar.extant]

current.data$All$info$Familia <- factor (current.data$All$info$Familia, levels = unique(current.data$All$info$Familia) )
current.data$All$info$Genero <- factor (current.data$All$info$Genero, levels = unique(current.data$All$info$Genero) )

contMap()
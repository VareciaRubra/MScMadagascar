####### Selecionando o que que vc quer plotar #########
current.data <- extant.madagascar.main.data

lin_data = ldply(current.data, function(x) tbl_df(cbind(x$info, x$ed))) 
lin_data %<>% mutate (., .info = paste(Familia, Genero, Especie, sep = "."))
individual_ID <- unlist(llply(current.data, function(x) x$info$Tombo))
myformula = paste0('~', paste(names(select(lin_data, IS_PM:BA_OPI)), collapse = '+'))
rownames(lin_data) <- individual_ID
PRCOMP <- princomp(as.formula(myformula), data = lin_data, scale = T)
PRCOMP$scores
######## Biplot PC1 x PC2 ##############
current.data <- select(lin_data, .info, IS_PM:BA_OPI)
rownames(current.data) <- rownames(lin_data)
#PRCOMP %>% biplot
#resp %<>% mutate(., PC1 = PRCOMP$scores[,1], PC2=PRCOMP$scores[,2]) 
ancestral.mx <- PhyloW(tree = pruned.tree.with.mx, tip.data = mx.at.tree, tip.sample.size = n.size[mask,2][-41])
plot(pruned.tree.with.mx, cex = 0.5)
nodelabels()
W.matrix <- ancestral.mx$'42'

resp <- as.data.frame(as.matrix(na.omit(select(current.data, IS_PM:BA_OPI) ) ) %*% eigen(W.matrix)$vectors[,c(1:39)])
resp %<>% mutate(., .id = rownames(resp) )
names(resp) %<>% gsub("V", "PC", .)
resp.info <- select(lin_data, c(.info,Arquivo:Data_dado) )
resp.info$Tombo <- reorder.factor(resp.info$Tombo, new.order=resp$.id)
resp.info <- resp.info %>% arrange(Tombo)
resp.info <- resp.info[resp.info$Tombo %in% resp$.id, ]
plot.W <- cbind(resp.info, resp)
plot.W$.id == plot.W$Tombo
rownames(plot.W) <- plot.W$Tombo
hulls <- ddply(plot.W, .(.info), plyr::summarise, "hpc1"=PC1[chull(PC1,PC2)],
               "hpc2"=PC2[chull(PC1,PC2)])
hulls %<>% separate(.info, c('Familia', 'Genero', "Especie"), sep = "\\.")

# levels(plot.W$Especie) <- as.factor(unique(resp.info$Especie))
# levels(plot.W$Genero) <- as.factor(unique(resp.info$Genero))
 levels(plot.W$Familia) <- as.factor(c("Indridae", "Daubentoniidae", "Lemuridae", "Lepilemuridae", "Cheirogaleidae", "Palaeopropithecidae", "Megaladapidae", "Archaeolemuridae", "Lorisidae", "Galagidae", "Tarsiidae" ))
# hulls$Genero <- factor (hulls$Genero, levels = hulls$Genero )
hulls$Especie <- unique (hulls$Especie, levels = unique(as.factor(hulls$Especie) ) )
pc_plot <- ggplot(plot.W, aes(PC1, PC2)) +
  geom_polygon(aes(-hpc1, -hpc2, fill = Familia, color = Familia, group= Genero), data = hulls, alpha=.3) + 
  geom_point(data = ddply(plot.W, .(Especie), numcolwise(mean)),
             aes(-PC1, -PC2, group= Especie), size = 3) + 
  theme_bw() + ggtitle("Cranial traits Within-group PC scores") + 
  coord_fixed()
pc_plot


hulls <- ddply(plot.W, .(.info), plyr::summarise, "hpc1"=PC1[chull(PC1,PC4)],
               "hpc2"=PC4[chull(PC1,PC4)])
hulls %<>% separate(.info, c('Familia', 'Genero', "Especie"), sep = "\\.")
hulls$Familia <- as.factor(hulls$Familia)
levels(hulls$Familia) <- order(unique(hulls$Familia)) as.factor(c("Indridae", "Daubentoniidae", "Lemuridae", "Lepilemuridae", "Cheirogaleidae"))
str(hulls)
pc_plot <- ggplot(plot.W, aes(PC1, PC4))+
  geom_polygon(aes(-hpc1, -hpc2, fill = Familia, color = Familia, group= Especie ), data = hulls, alpha=.2) + 
  geom_point(data = ddply(plot.W, .(Tombo), numcolwise(mean)),
             aes(-PC1, -PC4, group= Tombo), size = 1, color = "grey", alpha = 0.3) +
  #   geom_text(data = ddply(plot.W, .(Tombo), numcolwise(mean)),
  #              aes(PC1, PC2, label= Tombo), size = 3) 
  #   geom_point(data = ddply(plot.W, .(Genero), numcolwise(mean)),
  #             aes(PC1, PC2, group= Genero), size = 3) +
  geom_text(data = ddply(plot.W, .(Genero), numcolwise(mean)),
            aes(-PC1, -PC4, label= Genero), size = 5) +
  #scale_y_continuous(limits = c(-40, 15)) +
  #scale_x_continuous(limits = c(30, 190)) +
  theme_bw() + ggtitle("Cranial traits Within-group PC scores") +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 12),
        legend.position= c(0.37,0.17)) +
    guides(fill = guide_legend(keywidth = 1, keyheight = 2))
        #legend.position= "bottom") 
pc_plot


pc_plot <- ggplot(plot.W, aes(PC1, PC2)) +
  geom_polygon(aes(hpc1, hpc2, fill = Familia, color = Familia, group= Genero ), data = hulls, alpha=.2) + 
  geom_point(data = ddply(plot.W, .(Tombo), numcolwise(mean)),
             aes(PC1, PC2, group= Tombo), size = 1, color = "grey", alpha = 0.3) +
  #   geom_text(data = ddply(plot.W, .(Tombo), numcolwise(mean)),
  #              aes(PC1, PC2, label= Tombo), size = 3) 
  #   geom_point(data = ddply(plot.W, .(Genero), numcolwise(mean)),
  #             aes(PC1, PC2, group= Genero), size = 3) +
  geom_text(data = ddply(plot.W, .(Genero), numcolwise(mean)),
            aes(PC1, PC2, label= Genero), size = 5) +
  theme_bw() + ggtitle("Cranial traits Within-group PC scores")
pc_plot


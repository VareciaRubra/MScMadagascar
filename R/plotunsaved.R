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
W.matrix <- ancestral.mx$'45'

resp <- as.data.frame(as.matrix(na.omit(select(current.data, IS_PM:BA_OPI) ) ) %*% eigen(W.matrix)$vectors[,c(1:39)])
resp %<>% mutate(., .id = rownames(resp) )
names(resp) %<>% gsub("V", "PC", .)
resp.info <- select(lin_data, c(.info,Arquivo:Data_dado) )
resp.info$Tombo <- reorder.factor(resp.info$Tombo, new.order=resp$.id)
resp.info <- resp.info %>% arrange(Tombo)
resp.info <- resp.info[resp.info$Tombo %in% resp$.id, ]
table(resp$.id == resp.info$Tombo)
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
pc_plot + theme(legend.position = "none")


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

plot.W$PC1 <- -plot.W$PC1

plot.W.ranges <- plot.W
plot.W.ranges$Especie <- factor(plot.W.ranges$Especie, levels = unique(plot.W.ranges$Especie)[59:1] )

plot.W.ranges %>% gather(key=.pc, value = value, c(PC1, PC2, PC4, PC3 )) %>%
  ggplot( .,aes(y= Especie, x = value), varwidth = T) +
  geom_boxplot(aes(color = .pc, fill = .pc), alpha = 0.5, varwidth = T) +
  #geom_violin(aes(label = Tombo, color = .pc, fill = .pc), alpha = 0.5) +
  #geom_text( aes(label = Tombo), size =2, vjust = 1, alpha = 0.4)  +
  scale_colour_brewer(palette = "Spectral", direction = 1) +
  scale_fill_brewer(palette = "Spectral", direction = 1) +
  geom_jitter(aes(color = .pc), size = 0.3, shape = 8, alpha = 0.5) +
  #scale_shape( guide = "none", name = "Mean value by method") +
  ggtitle("Species scores on W's PCs") +
  theme_bw() + xlab("") + ylab("") +
  #coord_flip() +
  facet_wrap(~ .pc, ncol = 4, nrow = 1, scales = "free") +
  theme(legend.position="none",
        plot.title = element_text(lineheight=.8, face="bold"),
        axis.text.y = element_text(face =  "italic", size =13) 
        #axis.text.y = element_blank()
        )


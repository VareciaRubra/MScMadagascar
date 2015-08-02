
current.data <- extant.madagascar.main.data
current.data$All$info$Familia <- factor (current.data$All$info$Familia, levels = unique(current.data$All$info$Familia) )
current.data$All$info$Genero <- factor (current.data$All$info$Genero, levels = unique(current.data$All$info$Genero) )

lin_data = ldply(current.data, function(x) tbl_df(cbind(x$info, x$ed))) 
lin_data$Especie %<>% gsub("\\.", "", .)
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
Wmat = CalculateMatrix(manova(as.matrix(current.data$All$ed)  ~ Especie, data = as.data.frame(current.data$All$info) ) )
resp <- as.data.frame(as.matrix(na.omit(select(current.data$All$ed, IS_PM:BA_OPI) ) ) %*% eigen(Wmat)$vectors[,c(1:39)])
resp %<>% mutate(., .id = rownames(resp) )
names(resp) %<>% gsub("V", "PC", .)
resp.info <- select(lin_data, c(.info,Arquivo:Data_dado) )
resp.info <- resp.info[resp.info$Tombo %in% resp$.id, ]
plot.W <- cbind(resp.info, resp)
hulls <- ddply(plot.W, .(.info), plyr::summarise, "hpc1"=PC1[chull(PC1,PC4)],
               "hpc2"=PC4[chull(PC1,PC4)])
hulls %<>% separate(.info, c('Familia', 'Genero', 'Especie'), sep = "\\.")

pc_plot <- ggplot(plot.W, aes(PC1, PC4)) +
  geom_polygon(aes(hpc1, hpc2, fill = Familia, color = Familia, group= Genero), data = hulls, alpha=.3) + 
  geom_point(data = ddply(plot.W, .(Especie), numcolwise(mean)),
             aes(PC1, PC4, group= Especie), size = 3) + 
  theme_bw() + ggtitle("Cranial traits Within-group PC scores") + 
  coord_fixed()
pc_plot


pc_plot <- ggplot(plot.W, aes(PC1, PC4)) +
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
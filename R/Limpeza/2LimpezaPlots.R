############################  GRAFICOS PARA A LIMPEZA DE DADOS ########################
#1.Violin das medias das especies por variável
# raw.data %>%
#   #filter(., is.na(Sexo) ) %>% 
#   gather(key=ed, value=value, 12:50 ) %>% 
#   ggplot(., aes( x= Especie, y=value, color=Especie), varwidth = T) + 
#   geom_violin() + 
#   geom_jitter() +
#   facet_wrap(~ed, scale="free_y") + 
#   theme(axis.text.x = element_text(angle = 90)) 
  
# mesmo plot mas com os labels de cada indivíduo, bom para identificar os outliers e remove-los individualmente

raw.data %>%
  #filter(., Genero == "Eulemur" ) %>%
 #mutate(Loc.info = !is.na(Localidade)) %>%
  gather(key=ed, value=value, (16:54) ) %>% 
  ggplot(., aes( x= Museu , y=value, label = Tombo) ) + 
  geom_text(angle = 0, size =0.8)  +
  #geom_boxplot(aes(x= Especie , y=value), alpha = 0.5) + 
  geom_violin( aes(x= Museu, y=value), alpha = 0) + 
  geom_jitter(size = 0.8, aes(color = Museu), alpha = 0.3 )+
  #scale_shape(guide = "legend", name = "Specie") +
  facet_wrap(~ed,  scales="free") + 
  theme(axis.text.x = element_text(angle = 0)) +
  ggtitle("Traits distribuition") +
  theme(axis.title.x = element_blank() ) +
  #theme(legend= element_blank() ) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
  theme(plot.title = element_text(lineheight=.8, face="bold"))

###########################Vendo os caras que nao tem info de localidade
raw.data %>%
  filter(., Genero == "Microcebus" ) %>%
  mutate(Loc.info = !is.na(Localidade)) %>%
  gather(key=ed, value=value, (16:51) ) %>% 
  ggplot(., aes( x= interaction(Loc.info, Especie) , y=value, label = Tombo) ) + 
  geom_text(angle = 0, size =1)  +
  geom_boxplot(aes(x= interaction(Loc.info, Especie) , y=value, fill = Loc.info ), alpha = 0.5) + 
  #geom_violin( aes(x= Especie, y=value), alpha = 0) + 
  geom_jitter(aes(shape= Especie) )+
  scale_shape(guide = "legend", name = "Specie") +
  facet_wrap(~ed,  scales="free") + 
  #theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Species traits distribuition") +
  #theme(axis.title.x = element_blank() ) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank()) +
    theme(plot.title = element_text(lineheight=.8, face="bold"))
  

# mesmo plot com os labels de cada indivíduo mas por variável bom para identificar os outliers e remove-los individualmente
raw.data %>%
  #filter(., is.na(Sexo) ) %>% 
  gather(key=ed, value=value, IS_NSL ) %>% 
  ggplot(., aes( x= Especie, y=value, color=Especie, label = Tombo), varwidth = T) + 
  geom_text(aes(size = 0.4, vjust = 1) )  +
  geom_violin(aes(alpha = 0)) + 
  geom_jitter() +
  #facet_wrap(~ed, scale="free_y") + 
  theme(axis.text.x = element_text(angle = 90)) 


#ggplot da distribuiçao dos valires de determinada variável em determinada espécie num histograma
#ggplot(filter(raw.data, Especie == "Hapalemur_griseus"), aes(BA_OPI, color= Sexo)) + geom_histogram() + facet_grid(Especie~Sexo)

# #plot das medias de cada espécie por variável
# raw.data %>% 
#   dplyr::group_by(., Especie ) %>%
#   summarise_each(funs(var),IS_PM:BA_OPI) %>% 
#   gather(key=ed, value=value, 2:40 )%>%
#   ggplot(., aes( x= Especie, y=value, color=Especie)) +
#   geom_point(aes(group = Especie)) +
#   facet_wrap(~ed, scale="free_y") + 
#   theme(axis.text.x = element_text(angle = 90))

#plot das medias de cada variável por especie e ligando com um alinha 
#(pra ter uma noçao se tem algum a especie que é sempre maior/menor que as demais...)

raw.data %>% 
  dplyr::group_by(., Especie, Sexo ) %>%
  summarise_each(funs(mean),IS_PM:BA_OPI) %>% RColorBrewer
  gather(key=ed, value=value, 3:41 )  %>%
  ggplot(., aes( x= ed, y=value, color=Especie)) +
  geom_line(aes(group = Especie)) +
  facet_wrap(~Sexo)+
  theme(axis.text.x = element_text(angle = 90))

#plot em linhas por media entre replica de indivíduos por museu e por especie

raw.data %>% 
    dplyr::group_by(., Especie, Museu, Tombo ) %>%
  summarise_each(funs(mean),IS_PM:BA_OPI) %>% 
  gather(key=ed, value=value, 4:42 )  %>%
  ggplot(., aes( x= ed, y=value, color=Museu )) +
  geom_line(aes(group = Tombo, linetype = Especie)) +
  #facet_wrap(~Especie)+
  theme(axis.text.x = element_text(angle = 90))

##### com todo mundo #####
raw.data %>% 
  dplyr::group_by(., Familia, Genero, Especie, Museu, Tombo ) %>%
  summarise_each(funs(mean),IS_PM:BA_OPI) %>% 
  gather(key=ed, value=value, 6:42 )  %>%
  ggplot(., aes( x= ed, y=value, color=Genero )) +
  geom_line(aes(group = Tombo, linetype = Museu)) +
  facet_wrap(~Familia, scale="free_y")+
  theme(axis.text.x = element_text(angle = 90))

raw.data %>% 
  filter(Especie == "Propithecus_verreauxi") %>%
  dplyr::group_by(., Museu, Tombo ) %>%
  summarise_each(funs(mean),IS_PM:BA_OPI) %>% 
  gather(key=ed, value=value, 3:41 )  %>%
  ggplot(., aes( x= ed, y= value, color= Museu )) +
  geom_line(aes(group = Tombo)) +
  theme(axis.text.x = element_text(angle = 90))

######## Biplot PC1 x PC2 ##############
current.data <-extant.main.data$All
#PRCOMP <- princomp(current.data$ed %>% na.omit) # extraindo os componentes principais a partir dos dados (medias ed)
PRCOMP <- princomp(data.frame(na.omit(current.data$ed ) ) )
resp <- current.data$info[current.data$info$Tombo %in% dimnames(PRCOMP$scores)[[1]], ] #respectivos dados aos participantes da PCA
resp %<>% mutate(., PC1 = PRCOMP$scores[,1], PC2=PRCOMP$scores[,2]) 
hulls <-ddply(resp, .(Familia), plyr::summarise, "hpc1"=PC1[chull(PC1,PC2)],
              "hpc2"=PC2[chull(PC1,PC2)])
fig <- ggplot(resp, aes(PC1, PC2, color = Familia)) +
  geom_text(aes(PC1, PC2, color = Familia, label = Especie) ) +
  geom_polygon(aes(hpc1, hpc2, fill = Familia, group = Familia),
               hulls,alpha=.2)
fig


lin_data = ldply(extant.main.data, function(x) tbl_df(cbind(x$info, x$ed))) 
lin_data$Especie %<>% gsub("\\.", "", .)
lin_data %<>% mutate (., .info = paste(Familia, Genero, sep = "."))
individual_ID <- unlist(llply(extant.main.data, function(x) x$info$Tombo))
myformula = paste0('~', paste(names(select(lin_data, IS_PM:BA_OPI)), collapse = '+'))
rownames(lin_data) <- individual_ID
PRCOMP <- princomp(as.formula(myformula), data = lin_data, scale = T)
#PRCOMP$scores
######## Biplot PC1 x PC2 ##############
current.data <- select(lin_data, .info, IS_PM:BA_OPI)
rownames(current.data) <- rownames(lin_data)
#PRCOMP %>% biplot
#resp %<>% mutate(., PC1 = PRCOMP$scores[,1], PC2=PRCOMP$scores[,2]) 
Wmat = CalculateMatrix(manova(as.matrix(extant.main.data$All$ed)  ~ Especie, data = as.data.frame(extant.main.data$All$info) ) )
resp <- as.data.frame(as.matrix(na.omit(select(current.data, IS_PM:BA_OPI) ) ) %*% eigen(Wmat)$vectors[,1:2])
resp %<>% mutate(., ID = rownames(resp) )
names(resp) <- c("PC1", "PC2", ".id")
resp.info <- select(lin_data, c(.info,Arquivo:Data_dado) )
resp.info <- resp.info[resp.info$Tombo %in% resp$.id, ]
plot.W <- cbind(resp.info, resp)
hulls <- ddply(plot.W, .(.info), plyr::summarise, "hpc1"=PC1[chull(PC1,PC2)],
              "hpc2"=PC2[chull(PC1,PC2)])
hulls %<>% separate(.info, c('Familia', 'Genero'), sep = "\\.")

pc_plot <- ggplot(plot.W, aes(PC1, PC2)) +
  geom_polygon(aes(hpc1, hpc2, fill = Familia, color = Familia, group= Genero ), data = hulls, alpha=.3) + 
  geom_point(data = ddply(plot.W, .(Genero), numcolwise(mean)),
             aes(PC1, PC2, group= Genero), size = 3) + 
  geom_text(data = ddply(plot.W, .(Especie), numcolwise(mean)),
            aes(PC1, PC2, label= Especie), size = 3) +
  geom_point(data = ddply(plot.W, .(Tombo), numcolwise(mean)),
             aes(PC1, PC2, group= Tombo), size = 1) + 
             theme_bw() + ggtitle("Cranial traits Within-group PC scores")
pc_plot


pc_plot <- ggplot(plot.W, aes(PC1, PC2)) +
  geom_polygon(aes(hpc1, hpc2, fill = Familia, color = Familia, group= Especie ), data = hulls, alpha=.3) + 
  geom_text(data = ddply(plot.W, .(Especie), numcolwise(mean)),
            aes(PC1, PC2, label= Especie), size = 3) +
  theme_bw() + ggtitle("Cranial traits Within-group PC scores")
pc_plot

#################################### Sizeless ##################################

extant.madagascar.main.data$All$info$Familia <- factor (extant.madagascar.main.data$All$info$Familia, levels = unique(extant.madagascar.main.data$All$info$Familia) )
extant.madagascar.main.data$All$info$Genero <- factor (extant.madagascar.main.data$All$info$Genero, levels = unique(extant.madagascar.main.data$All$info$Genero) )

lin_data = ldply(extant.madagascar.main.data, function(x) tbl_df(cbind(x$info, x$ed))) 
lin_data$Especie %<>% gsub("\\.", "", .)
lin_data %<>% mutate (., .info = paste(Familia, Genero, sep = "."))
individual_ID <- unlist(llply(extant.madagascar.main.data, function(x) x$info$Tombo))
myformula = paste0('~', paste(names(select(lin_data, IS_PM:BA_OPI)), collapse = '+'))
rownames(lin_data) <- individual_ID
PRCOMP <- princomp(as.formula(myformula), data = lin_data, scale = T)
#PRCOMP$scores
######## Biplot PC1 x PC2 ##############
current.data <- select(lin_data, .info, IS_PM:BA_OPI)
rownames(current.data) <- rownames(lin_data)
#PRCOMP %>% biplot
#resp %<>% mutate(., PC1 = PRCOMP$scores[,1], PC2=PRCOMP$scores[,2]) 
Wmat = CalculateMatrix(manova(as.matrix(extant.madagascar.main.data$All$ed)  ~ Especie, data = as.data.frame(extant.madagascar.main.data$All$info) ) )
resp <- as.data.frame(as.matrix(na.omit(select(current.data, IS_PM:BA_OPI) ) ) %*% eigen(Wmat)$vectors[,c(1:39)])
resp %<>% mutate(., .id = rownames(resp) )
names(resp) %<>% gsub("V", "PC", .)
resp.info <- select(lin_data, c(.info,Arquivo:Data_dado) )
resp.info <- resp.info[resp.info$Tombo %in% resp$.id, ]
plot.W <- cbind(resp.info, resp)
hulls <- ddply(plot.W, .(.info), plyr::summarise, "hpc1"=PC1[chull(PC1,PC2)],
               "hpc2"=PC2[chull(PC1,PC2)])
hulls %<>% separate(.info, c('Familia', 'Genero'), sep = "\\.")

pc_plot <- ggplot(plot.W, aes(PC1, PC2)) +
  geom_polygon(aes(hpc1, hpc2, fill = Familia, color = Familia, group= Genero), data = hulls, alpha=.3) + 
  geom_point(data = ddply(plot.W, .(Especie), numcolwise(mean)),
             aes(PC1, PC2, group= Especie), size = 3) + 
  theme_bw() + ggtitle("Cranial traits Within-group PC scores") + 
  coord_fixed()
pc_plot


pc_plot <- ggplot(plot.W, aes(PC1, PC2)) +
  geom_polygon(aes(hpc1, hpc2, fill = Familia, color = Familia, group= Especie ), data = hulls, alpha=.2) + 
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



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
  #filter(., is.na(Sexo) ) %>% 
  gather(key=ed, value=value, 12:50 ) %>% 
  ggplot(., aes( x= Especie, y=value, color=Especie, label = Tombo), varwidth = T) + 
  geom_text(aes(size =.4, vjust = 1) )  +
  geom_violin(aes(alpha = 0)) + 
  geom_jitter() +
  facet_wrap(~ed, scale="free_y") + 
  theme(axis.text.x = element_text(angle = 90)) 

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
  summarise_each(funs(mean),IS_PM:BA_OPI) %>% 
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

raw.data %>% 
  filter(Especie == "Propithecus_verreauxi") %>%
  dplyr::group_by(., Museu, Tombo ) %>%
  summarise_each(funs(mean),IS_PM:BA_OPI) %>% 
  gather(key=ed, value=value, 3:41 )  %>%
  ggplot(., aes( x= ed, y= value, color= Museu )) +
  geom_line(aes(group = Tombo)) +
  theme(axis.text.x = element_text(angle = 90))

######## Biplot PC1 x PC2 ##############
current.data <- all.main.data$Lemur
PRCOMP <- princomp(current.data$ed) # extraindo os componentes principais a partir dos dados (medias ed)
#PRCOMP <- princomp(data.frame(na.omit(current.data$sizeless ) ) )
resp <- current.data$info[current.data$info$Tombo %in% dimnames(PRCOMP$scores)[[1]], ] #respectivos dados aos participantes da PCA
resp %<>% mutate(., PC1 = PRCOMP$scores[,1], PC2=PRCOMP$scores[,2]) 
hulls <-ddply(resp, .(Especie), plyr::summarise, "hpc1"=PC1[chull(PC1,PC2)],
              "hpc2"=PC2[chull(PC1,PC2)])
fig <- ggplot(resp, aes(PC1, PC2, color= Especie)) +
  geom_text(aes(PC1, PC2, color=Especie, label = Tombo) ) +
  geom_polygon(aes(hpc1, hpc2, fill= Especie, group= Especie),
               hulls,alpha=.2)
fig




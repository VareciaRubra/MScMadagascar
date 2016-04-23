Modular.hyp <- Aux$def.hyp[c(2:20, 1, 21:39),] # pegando as hipoteses modulares do arquivo auxiliat do Pato
Modular.hyp[20,] <- c(0,0,1,0,0,1,1,1) # mas ele tem logCS e nao tem PT.TSP, enfiemolhos!
rownames(Modular.hyp)[20] <- "PT.TSP"

Modulemurs <- vector("list")
Modulemurs$test.modularity <- sp.main.data[mask.no.na.cov] %>% 
  llply (function (x) TestModularity(cor.matrix = x$matrix$cor, modularity.hipot = Modular.hyp , permutations = 1000, MHI = F), .progress = "text")

Modulemurs$Probability<- Modulemurs$test.modularity %>% ldply(function (x) x$Probability)
names(Modulemurs$Probability) <- c("sp", Modulemurs$test.modularity$Tarsius_bancanus$hypothesis)
Modulemurs$Probability$sp <- factor (Modulemurs$Probability$sp, levels = rev(unique(Modulemurs$Probability$sp)))
Modulemurs$Probability$sp %<>% gsub("_", " ",.)
Modulemurs$Probability$sp <- factor (Modulemurs$Probability$sp, levels = rev(unique(Modulemurs$Probability$sp)))

Modulemurs$Probs<- Modulemurs$test.modularity %>% ldply(function (x) x$Probability)
names(Modulemurs$Probs) <- c("sp", Modulemurs$test.modularity$Tarsius_bancanus$hypothesis)
rownames(Modulemurs$Probs) <- Modulemurs$Probs[,1]
Modulemurs$Probs[Modulemurs$Probs <= 0.05] <- "red"
Modulemurs$Probs[Modulemurs$Probs != "red"] <- "darkgrey"
Modulemurs$Probs[,1] <- rownames(Modulemurs$Probs)

Modulemurs$Probability %>% gather(key = "Probability", value = probability, 2:10) %>%
ggplot (.) +
  geom_tile(aes(x = Probability, y = sp, fill = as.numeric(probability<= 0.05) ), alpha = 0.9,  color = "grey") +
  #geom_point(aes(x =Probability, y = sp, color = as.numeric(probability<= 0.05) ) ) +
  theme_bw() +
  scale_fill_gradient(name = 'Probability', low = "transparent", high = "red") +
  #scale_color_gradient(name = 'Probability', low = "transparent", high = "red" ) +
  ylab ('') + xlab ('') +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1, size = 8),
        axis.ticks = element_line(size = 0),
        axis.text.y = element_text(face = "italic", size= 10),
        legend.position = "none",
        rect = element_blank(), 
        line = element_blank())

Modulemurs$AVG.ratio<- Modulemurs$test.modularity %>% ldply(function (x) x$`AVG Ratio`)
names(Modulemurs$AVG.ratio) <- c("sp", Modulemurs$test.modularity$Tarsius_bancanus$hypothesis)
Modulemurs$AVG.ratio$sp <- factor (Modulemurs$AVG.ratio$sp, levels = rev(unique(Modulemurs$AVG.ratio$sp)))
Modulemurs$AVG.ratio$sp %<>% gsub("_", " ",.)
Modulemurs$AVG.ratio$sp <- factor (Modulemurs$AVG.ratio$sp, levels = rev(unique(Modulemurs$AVG.ratio$sp)))
rownames(Modulemurs$AVG.ratio) <- Modulemurs$AVG.ratio[,1]
Modulemurs$AVG.ratio[,1] <- rownames(Modulemurs$AVG.ratio)

Modulemurs$AVG.ratio %>% gather(key = "Region", value = AVG.ratio, 2:10) %>% 
  ggplot (.) +
  geom_tile(aes(x = Region, y = sp, fill = abs(AVG.ratio)) , alpha = 0.9,  color = "grey") +
  #geom_point(aes(x =AVG.ratio, y = sp, color = as.numeric(AVG.ratio<= 0.05) ) ) +
  theme_bw() +
  scale_fill_gradient2(name = 'AVG ratio', low = "transparent", mid = "#f46d43", high ="#d73027", midpoint = 3) +
  ylab ('') + xlab ('') +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1, size = 8),
        axis.ticks = element_line(size = 0),
        axis.text.y = element_text(face = "italic", size= 10),
        #legend.position = "none",
        rect = element_blank(), 
        line = element_blank())
PC1.percent$.id %in% criminosos
criminosos_ <- c("Euoticus_elegantulus", "Prolemur_simus", "Tarsius_syrichta")
criminosos = PC1.percent$.id %in% criminosos_
criminosos <- criminosos== FALSE
Modulemurs$cor.pc1.flex <- cbind.data.frame("cor.Flex" = cor(abs(Modulemurs$AVG.ratio[criminosos, 2:10]), Flex$V1[criminosos]), 
                                            "cor.PC1" = cor(abs(Modulemurs$AVG.ratio[criminosos, 2:10]), PC1.percent$V1[criminosos]),
                                            "cor.r2" = cor(abs(Modulemurs$AVG.ratio[criminosos, 2:10]), r2.stuff[42:1,2][criminosos,]))
r2.stuff[42:1,]
Flex <- Variae %>% ldply( function(x) mean(x$flex.dist))
PC1.percent <- Variae %>% ldply( function(x) x$intervalo.mc.pc$observed[1] )
R2 <- Variae %>% ldply(function(x) mean(x$r2.dist ))
AVG.r <- Modulemurs$AVG.ratio %>% gather(key = "Region", value = AVG.ratio, 2:10) 
AVG.r$Flex <- rep(Flex$V1, 9) 
AVG.r$PC1.percent <- rep(PC1.percent$V1, 9) 
AVG.r$r2 <- rep(R2$V1, 9) 
AVG.r$cor.flex <- rep(Modulemurs$cor.pc1.flex$cor.Flex, each = length(unique(AVG.r$sp)) )
AVG.r$cor.PC1<- rep(Modulemurs$cor.pc1.flex$cor.PC1, each = length(unique(AVG.r$sp)) )
AVG.r$cor.r2<- rep(Modulemurs$cor.pc1.flex$mean, each = length(unique(AVG.r$sp)) )


Modulemurs$Plot$AVG.Flex <- AVG.r %>% tbl_df %>% filter (sp != "Euoticus elegantulus") %>% filter (sp != "Prolemur simus")%>% filter (sp != "Tarsius syrichta")%>%
  ggplot (aes(group = Region)) +
 geom_point(aes(y = abs(AVG.ratio), x = Flex), size = 1, color = "darkgrey") + 
  geom_text(aes(y = 4.2, x = 0.37 , label = paste ("r=", round(cor.flex,2) ), group = Region)) +
  geom_text(aes(y = AVG.ratio, x = Flex, label = sp), size = 3, alpha = 0.8) + 
  geom_smooth(aes(y = abs(AVG.ratio), x = Flex), method = "lm") +
  facet_wrap(~Region) + theme_bw() +
  theme(axis.text = element_text(size = 7))
Modulemurs$Plot$AVG.PC1 <- AVG.r%>% tbl_df %>% filter (sp != "Euoticus elegantulus") %>% filter (sp != "Prolemur simus")%>% filter (sp != "Tarsius syrichta")%>%
  ggplot (aes(group = Region )) +
  geom_point(aes(y = abs(AVG.ratio), x = PC1.percent), size = 1, color = "darkgrey") + 
  geom_text(aes(y = 4.2, x = 0.45 , label = paste("r=", round(cor.PC1, 2)), group = Region)) +
  geom_text(aes(y = AVG.ratio, x = PC1.percent, label = sp),  size = 3, alpha = 0.8) + 
  geom_smooth(aes(y = abs(AVG.ratio), x = PC1.percent), method = "lm", alpha = 0.5) +
  facet_wrap(~Region)+ theme_bw() +
  theme(axis.text = element_text(size = 7))

plot_grid(Modulemurs$Plot$AVG.Flex, Modulemurs$Plot$AVG.PC1, ncol = 2) #########################################################################################

Modulemurs$Plot$AVG.R2 <- AVG.r %>% tbl_df %>% filter (sp != "Euoticus elegantulus") %>% filter (sp != "Prolemur simus")%>% filter (sp != "Tarsius syrichta")%>%
  ggplot (aes(group = Region)) +
  geom_point(aes(y = abs(AVG.ratio), x = r2), size = 1, color = "darkgrey") + 
  geom_text(aes(y = 4.2, x = 0.22 , label = paste ("r=", round(cor.r2,2) ), group = Region)) +
  geom_text(aes(y = AVG.ratio, x = r2, label = sp), size = 3, alpha = 0.8) + 
  geom_smooth(aes(y = abs(AVG.ratio), x = r2), method = "lm") +
  facet_wrap(~Region) + theme_bw() +
  theme(axis.text = element_text(size = 7))

DistModular <- function (x, simulations = 1000, modularity.hipot = Modular.hyp)
  {
  hypothesis <- c(dimnames(Modular.hyp)[[2]], "Full.integration")
  Rsquared <- matrix(data = NA, nrow = simulations, ncol = dim(Modular.hyp)[[2]] +1)
  colnames(Rsquared) <- hypothesis
  AVG.plus <- matrix(data = NA, nrow = simulations, ncol = dim(Modular.hyp)[[2]] + 1)
  colnames(AVG.plus) <- hypothesis
  AVG.minus <- matrix(data = NA, nrow = simulations, ncol = dim(Modular.hyp)[[2]]+ 1)
  colnames(AVG.minus) <- hypothesis
  Probability <- matrix(data = NA, nrow = simulations, ncol = dim(Modular.hyp)[[2]]+ 1)
  colnames(Probability) <- hypothesis
  AVG.rat <- matrix(data = NA, nrow = simulations, ncol = dim(Modular.hyp)[[2]]+ 1)
  colnames(AVG.rat) <- hypothesis
 
  observed <- TestModularity(cor.matrix = x$matrix$cor, modularity.hipot = Modular.hyp , permutations = 1, MHI = F)
  
  for (i in 1:simulations)  {
    Mx.Cor <- var (rmvnorm (x$sample.size, sigma = x$matrix$cor, method = 'svd') )
    modular <- TestModularity(cor.matrix = Mx.Cor, modularity.hipot = Modular.hyp , permutations = 1, MHI = F) 
    
    Rsquared [i, ] <- modular$Rsquared
    AVG.plus [i, ] <- modular$`AVG+`
    AVG.minus [i, ] <- modular$`AVG-`
    Probability [i, ] <-  modular$Probability
    AVG.rat [i, ] <-  modular$`AVG Ratio`
  }
  plotaisso <- cbind(rep (unique(x$info.raw$Especie), simulations),
  #melt(Probability, value.name = "Probability")[2:3],
  melt(AVG.plus, value.name = "AVG+")[2:3],
  melt(AVG.minus, value.name = "AVG-")[3] )
  names(plotaisso)[1:2] <- c("sp", "hyp")
  
  who <- unique(x$info.raw$Especie)
  #Probs <- cbind(rep (unique(x$info.raw$Especie), simulations), melt(Probability, value.name = "Probability")[2:3])
  #names(Probs)[1:2] <- c("sp", "hyp")
  AVG.r <- cbind(rep (unique(x$info.raw$Especie), simulations), melt(AVG.rat, value.name = "AVG.rat")[2:3])
  names(AVG.r)[1:2] <- c("sp", "hyp")
  # AVG.r.plot <- AVG.r %>% 
  #   ggplot(., aes (x = hyp, y = abs(AVG.rat), group = hyp)) +
  #   geom_violin() +
  #   ggtitle(who) 

  Plotavg <- data.frame("sp" = as.character(rep(unique(AVG.r$sp), 9)),
        "hypotesis"= observed$hypothesis,
        "AVG.ratio" = observed$`AVG Ratio`,
        "mean" = AVG.r %>%  group_by(hyp) %>% summarise_each(funs(mean(abs(.))), abs(AVG.rat)) %>% as.data.frame%>% .[,2],
        "min" = AVG.r %>%  group_by(hyp) %>% summarise_each(funs(sort(abs(.))[simulations*0.025]), abs(AVG.rat) ) %>% as.data.frame%>% .[,2],
        "max" = AVG.r %>%  group_by(hyp) %>% summarise_each(funs(sort(abs(.))[simulations*0.975]), abs(AVG.rat) ) %>% as.data.frame %>% .[,2] )
  
  AVG.r.plot <- Plotavg %>% 
    ggplot (data = ., aes (group = sp) ) + 
    geom_linerange(aes(x = sp, ymin = min, ymax = max, group = interaction (hypotesis, sp)), size =4, alpha = 0.2) +
    geom_point(aes(x = sp, y = mean, group = interaction (hypotesis, sp)), size = 4, color = "darkgrey", alpha = 0.7) +
    geom_point(aes(x = sp, y = AVG.ratio, group = interaction (hypotesis, sp))) +
    facet_wrap(~ hypotesis) + theme_bw() + coord_flip() + xlab("") + ylab("AVG ratio")

plot <- plotaisso %>% melt %>%
    ggplot(., aes(x = variable, y = abs(value), group = interaction(hyp, variable))) +
    geom_boxplot(aes(group = variable)) +
    ggtitle(who) +
    facet_wrap(~hyp, ncol = 3) + theme_bw() + ylab("absolute value") + xlab("")

plotaisso <- cbind(rep (unique(x$info.raw$Especie), simulations),
                   melt(Probability, value.name = "Probability")[2:3],
                   melt(AVG.plus, value.name = "AVG+")[2:3],
                   melt(AVG.minus, value.name = "AVG-")[3] )
names(plotaisso)[1:2] <- c("sp", "hyp")

  return(list ("plotaisso" = Plotavg,
               "AVG.r" =  AVG.r,
               "AVG.+-.plot" = plot,
               "AVG.r.plot" = AVG.r.plot
                ) )
}

temp <- DistModular(x = sp.main.data$Varecia_variegata, simulations = 10, modularity.hipot =  Modular.hyp )
temp$AVG.r.plot
temp$`AVG.+-.plot`
temp$AVG.r %>% ggplot(., aes (x = hyp, y = AVG.rat, group = hyp)) +
  geom_violin() +
  geom_boxplot() +
  geom_point(aes (x = hyp, y = mean(AVG.rat), group = hyp) ) +
  #geom_errorbar(aes(x = hyp, ymin = min(AVG.r), ymax = max(AVG.r), group = hyp ))+
  ggtitle(who) 

Modulemurs$test.modularity.dist <- sp.main.data[mask.no.na.cov] %>% llply(function (x) DistModular(x = x, simulations = 1000, modularity.hipot =  Modular.hyp) , .progress = "text")
save(Modulemurs, file = "Data/Modularity.RData")

Modulemurs$Plot$Oral <- Modulemurs$test.modularity.dist %>%  # pot em 10 x 15
  ldply(function(x) x$plotaisso) %>% filter (sp != "Euoticus_elegantulus") %>% 
  filter (sp != "Prolemur_simus") %>% 
  filter (hypotesis == "Oral") %>% 
  ggplot (data = ., aes (group = sp) ) + 
  geom_linerange(aes(x = sp, ymin = min, ymax = max, group = interaction (hypotesis, sp)), size =4, alpha = 0.2) +
  geom_point(aes(x = sp, y = mean, group = interaction (hypotesis, sp)), size = 4, color = "darkgrey", alpha = 0.7) +
  geom_point(aes(x = sp, y = abs(AVG.ratio), group = interaction (hypotesis, sp))) +
  facet_wrap(~ hypotesis) + theme_bw() + coord_flip() + xlab("") + ylab("AVG ratio")

Modulemurs$Plot$AVG.ratio <- plot_grid(Modulemurs$Plot$Oral, 
          Modulemurs$Plot$Nasal + theme(axis.text.y = element_blank() ), 
          Modulemurs$Plot$Zygomatic+ theme(axis.text.y = element_blank() ), 
          Modulemurs$Plot$Orbit, 
          Modulemurs$Plot$Base+ theme(axis.text.y = element_blank() ), 
          Modulemurs$Plot$Vault+ theme(axis.text.y = element_blank() ), 
          Modulemurs$Plot$Face, 
          Modulemurs$Plot$Neuro+ theme(axis.text.y = element_blank() ), 
          Modulemurs$Plot$Full+ theme(axis.text.y = element_blank() ), ncol = 3, rel_widths = c(1.8,1,1,1.8,1,1,1.9,1,1))  ################ 

Modulemurs$test.modularity.dist %>% ldply(function(x) x$plotaisso) %>% melt %>% 
  ggplot(., aes(x = variable, y = value, group = interaction(variable, hypotesis, .id))) +
  geom_boxplot(aes(group = interaction (variable, hypotesis, .id))) +
  facet_wrap(~hypotesis, ncol = 3)

Modulemurs$test.modularity.dist$Hapalemur_griseus$`AVG.+-.plot`

Modulemurs$Plot$AVG.ratio.selected <- Modulemurs$test.modularity.dist %>%  # pot em 10 x 15
  ldply(function(x) x$plotaisso) %>%  tbl_df %>%
  filter (sp == "Galago_senegalensis" | 
            sp == "Otolemur_crassicaudatus" |
            sp == "Euoticus_elegantulus" |
            sp == "Nycticebus_coucang" |
            sp == "Loris_tardigradus" | 
            sp == "Perodicticus_potto" |
            sp == "Daubentonia_madagascariensis" |
            sp == "Varecia_variegata" | 
            sp == "Lemur_catta" | 
            sp == "Hapalemur_griseus" | 
            sp == "Eulemur_rubriventer" |
            sp == "Indri_indri" | 
            sp == "Propithecus_verreauxi" | 
            sp == "Avahi_laniger" | 
            sp == "Lepilemur_mustelinus" | 
            sp == "Phaner_furcifer" |
            sp == "Cheirogaleus_medius" | 
            sp == "Mirza_coquereli" | 
            sp == "Microcebus_griseorufus" |
            sp == "Tarsius_bancanus")

Modulemurs$Plot$AVG.ratio.selected$hypotesis <- factor (Modulemurs$Plot$AVG.ratio.selected$hypotesis, levels = c(dimnames(Modular.hyp)[[2]], "Full Integration") ) 
Modulemurs$Plot$AVG.ratio.selected$sp <- factor (Modulemurs$Plot$AVG.ratio.selected$sp, levels = unique(Modulemurs$Plot$AVG.ratio.selected$sp) ) 
Modulemurs$Plot$AVG.ratio.selected$sp %<>% gsub("_", " ", .)
Modulemurs$Plot$AVG.ratio.selected$sp <- factor (Modulemurs$Plot$AVG.ratio.selected$sp, levels = unique(Modulemurs$Plot$AVG.ratio.selected$sp) ) 

Modulemurs$Plot$AVG.ratio.selected %>% 
  ggplot (data = ., aes (group = sp) ) + 
  geom_linerange(aes(x = sp, ymin = min, ymax = max, group = interaction (sp, hypotesis)), size =4, alpha = 0.2) +
  geom_point(aes(x = sp, y = mean, group = interaction (sp, hypotesis)), size = 4, color = "darkgrey", alpha = 0.7) +
  geom_point(aes(x = sp, y = abs(AVG.ratio), group = interaction (sp, hypotesis))) +
  facet_wrap(~ hypotesis) + theme_bw() + coord_flip() + xlab("") + ylab("AVG ratio") + 
  theme(axis.text = element_text(size = 17), 
        strip.text =element_text(size = 21),
        axis.text.y = element_text(face = "italic"))

Modulemurs$Plot$AVG.mm.selected <- Modulemurs$test.modularity.dist %>%  # pot em 10 x 15
  llply(function(x) x$`AVG.+-.plot`)

plot_grid(Modulemurs$Plot$AVG.mm.selected$Galago_senegalensis,
          Modulemurs$Plot$AVG.mm.selected$Otolemur_crassicaudatus,
          Modulemurs$Plot$AVG.mm.selected$Euoticus_elegantulus,
          Modulemurs$Plot$AVG.mm.selected$Nycticebus_coucang,
          Modulemurs$Plot$AVG.mm.selected$Loris_tardigradus,
          Modulemurs$Plot$AVG.mm.selected$Perodicticus_potto,
          Modulemurs$Plot$AVG.mm.selected$Daubentonia_madagascariensis,
          Modulemurs$Plot$AVG.mm.selected$Varecia_variegata,
          Modulemurs$Plot$AVG.mm.selected$Lemur_catta,
          Modulemurs$Plot$AVG.mm.selected$Hapalemur_griseus,
          Modulemurs$Plot$AVG.mm.selected$Eulemur_rubriventer,
          Modulemurs$Plot$AVG.mm.selected$Indri_indri,
          Modulemurs$Plot$AVG.mm.selected$Propithecus_verreauxi,
          Modulemurs$Plot$AVG.mm.selected$Avahi_laniger,
          Modulemurs$Plot$AVG.mm.selected$Lepilemur_mustelinus,
          Modulemurs$Plot$AVG.mm.selected$Phaner_furcifer,
          Modulemurs$Plot$AVG.mm.selected$Cheirogaleus_medius,
          Modulemurs$Plot$AVG.mm.selected$Mirza_coquereli,
          Modulemurs$Plot$AVG.mm.selected$Microcebus_griseorufus,
          Modulemurs$Plot$AVG.mm.selected$Tarsius_bancanus, ncol =3)

Modulemurs$Plot$AVG.ratio.selected <- Modulemurs$test.modularity.dist %>%  # pot em 10 x 15
  ldply(function(x) x$plotaisso) %>%  tbl_df %>%
  filter (sp == "Hapalemur_griseus"  )
            
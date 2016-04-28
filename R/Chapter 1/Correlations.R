

mean.sim<-function(x) {
  x[lower.tri(x)]<-t(x)[lower.tri(x)]
  diag(x)<-NA
  x<-rowMeans(x, na.rm=TRUE)
  return(x)
}

mean.sd<-function(x) {
  x[lower.tri(x)]<-t(x)[lower.tri(x)]
  diag(x)<-NA
  x<-apply(x, 1, sd, na.rm = TRUE)
  return(x)
}

mean.sim(abs(sp.main.data$Tarsius_bancanus$matrix$cor))

Moduling <- mx.list.taxonomy[-c(73,74)] %>% ldply( function(x) data.frame( "trait" =  names(mean.sim(x)), "mean.cor" = mean.sim(abs(cov2cor(x))), "sd.cor" = mean.sd(abs(cov2cor(x))) ) , .progress = "text" )

Moduling <- sp.main.data[mask.n.size] %>% 
  ldply( function(x) data.frame( "trait" =  names(mean.sim(x$matrix$cor)), "mean.cor" = mean.sim(abs(x$matrix$cor)), "sd.cor" = mean.sd(abs(x$matrix$cor)) ) , .progress = "text" )


Moduling$cv <- Moduling$sd.cor/ Moduling$mean.cor

Moduling$trait <- factor(Moduling$trait, levels = rev(unique(Moduling$trait)) )
Moduling$.id <- factor(Moduling$.id, levels = rev(unique(Moduling$.id)) )   
Moduling$r2 <- 
  rep(r2.stuff$mean, each = 39)

Moduling$gm <- 
  rep(gm.mean$V1, each = 39)
dimnames(Modular.hyp)
unique(Moduling$trait)

modules <- as.matrix(apply(Modular.hyp,2, FUN = rep, dim(Moduling)[1]/39))
dimnames(modules)
modules[,1]  
Moduling <- cbind(Moduling, modules)
Moduling$Neuro

names(Moduling)[1] <- "sp"
Moduling %>% filter (sp != "Euoticus_elegantulus") %>% 
  filter (sp != "Prolemur_simus") %>%
ggplot( aes(x = mean.cor, y = sd.cor)) +
  geom_point() +
  geom_text(aes(label = sp))+
  geom_smooth(method = "lm", se = F, aes(x = mean.cor, y = sd.cor, color = interaction (Face, Neuro), group = interaction (trait))) +
  facet_wrap(~trait, ncol = 8) +
  theme(legend.position = "none")

Moduling %>% 
  ggplot(.) +
  geom_density(aes(x= cv, group = trait)) +
  geom_histogram(aes(x= cv, group = trait, color = Base)  )  +
  facet_wrap(~trait)+ theme_bw()

Moduling %>% #filter (Face == 1) %>%
ggplot(., aes(y = cv, x = sp)) +
  geom_jitter(aes(color = trait))+
  geom_text(aes(label = trait)) +
   #geom_dotplot(stackdir = "up", stackratio = 0.3)+
  #facet_wrap(~trait) + 
  guides(color = guide_legend(ncol =1)) +theme(axis.text.x = element_text(angle = 90))


Moduling %>% #filter (Face == 1) %>%
  ggplot(., aes(y = cv, x = sp)) +
  geom_jitter()+
  geom_text(aes(label = trait, color = interaction(Vault, Neuro) ), size = 3) +
  #geom_dotplot(stackdir = "up", stackratio = 0.3)+
  #facet_wrap(~trait) + 
  guides(color = guide_legend(ncol =1)) +theme(axis.text.x = element_text(angle = 90))

Moduling %>% #filter (trait == "BA_OPI") %>%
  ggplot( aes(y = r2, x = cv)) +
  geom_point(aes(color = interaction (as.factor(Face), as.factor(Neuro)) )) +
  #geom_text(aes(label = trait, color = as.factor(Vault)) )+
  #geom_smooth(method = "lm", se = F, aes(x = cv, y = mean.cor, group = interaction (trait))) +
  theme_bw() + 
  facet_wrap(~trait, ncol = 8) 


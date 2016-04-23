

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

Moduling <- sp.main.data[mask.no.na.cov] %>% 
  ldply( function(x) data.frame( "trait" =  names(mean.sim(x$matrix$cor)), "mean.cor" = mean.sim(abs(x$matrix$cor)), "sd.cor" = mean.sd(abs(x$matrix$cor)) ) , .progress = "text" )

Moduling$cv <- Moduling$sd.cor/ Moduling$mean.cor

Moduling$trait <- factor(Moduling$trait, levels = rev(unique(Moduling$trait)) )
Moduling$.id <- factor(Moduling$.id, levels = rev(unique(Moduling$.id)) )   
Moduling$r2 <- 
  rep(r2.stuff$mean, each = 39)

modules <- as.matrix(apply(Modular.hyp,2, FUN = rep, 42))
dimnames(modules)
modules[,1]  
Moduling <- cbind(Moduling, modules)
Moduling$Neuro

names(Moduling)[1] <- "sp"
Moduling %>%   
ggplot( aes(x = mean.cor, y = sd.cor, color = trait)) +
  geom_smooth(method = "lm", se = F, aes(x = mean.cor, y = sd.cor, group = interaction (trait))) +
  facet_wrap(~trait, ncol = 8) +
  theme(legend.position = "none")

Moduling %>% melt %>% filter (variable == "cv") %>%
  ggplot( aes(y = trait, x = value)) +
  geom_violin() +
  geom_jitter(aes(color = as))
  geom_boxplot()+
  facet_wrap(~sp)

Moduling %>% #filter (trait == "BA_OPI") %>%
  ggplot( aes(y = r2, x = cv)) +
  geom_point(aes(color = interaction (as.factor(Base), as.factor(Vault)) )) +
  #geom_text(aes(label = trait, color = as.factor(Vault)) )+
  #geom_smooth(method = "lm", se = F, aes(x = cv, y = mean.cor, group = interaction (trait))) +
  theme_bw() + 
  facet_wrap(~trait, ncol = 8) + 


r2.stuff
Moduling %>% names


llply()
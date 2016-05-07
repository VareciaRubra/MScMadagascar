


SRD.Tree <- function (tree, cov.matrix.list, sample.sizes = NULL, mean.list) 
{
  if (!all(tree$tip.label %in% names(cov.matrix.list))) 
    stop("All tip labels must be in names(cov.matrix.list).")
  cov.matrices <- PhyloW(tree, cov.matrix.list, sample.sizes)
  
  nodes <- names(cov.matrices)
  node.mask <- laply(nodes, function(x) length(getMeans(mean.list, tree, x))) > 3
  
  newvec <- names(cov.matrices[1:70])
  for (i in unique(tree$tip.label)) newvec[names(cov.matrices[1:70]) == i ] <- c(1:70)[tree$tip.label == i ]
  names(cov.matrices)[1:70] <- newvec
  
  SRD.by.node <- llply(nodes[node.mask], 
                       function(node) 
                        {
                          indices <- as.character(tree$edge[,2] [tree$edge[,1] == node])
                          SRD.result <- tryCatch (expr = SRD(cov.matrices [[ indices [1] ]], cov.matrices [[ indices [2] ]]) , error = function(cond) vector("list") )
<<<<<<< HEAD
                          m.sd.cor <-list(abs(cov2cor(cov.matrices [[ indices [1] ]])), abs(cov2cor(cov.matrices [[ indices [2] ]])) ) %>% llply( function(x) data.frame( "trait" =  names(mean.sim(x)), "mean.cor" = mean.sim(abs(x)), "sd.cor" = mean.sd(abs(x)) ) , .progress = "text" )
                          RS.result <- RandomSkewers(cov.matrices [[ indices [1] ]], cov.matrices [[ indices [2] ]]) 
                          return(list(SRD.result = SRD.result, RS.result = RS.result, m.sd.cor = m.sd.cor))
=======
                          #m.sd.cor <-list(abs(cov2cor(cov.matrices [[ indices [1] ]])), abs(cov2cor(cov.matrices [[ indices [2] ]])) ) %>% llply( function(x) data.frame( "trait" =  names(mean.sim(x)), "mean.cor" = mean.sim(abs(x)), "sd.cor" = mean.sd(abs(x)) ) , .progress = "text" )
                          return(list(SRD.result = SRD.result))
>>>>>>> a6f0b3fa071439c4d680c4d364e9bebc4f3ed6d4
                              })
  names(SRD.by.node) <- paste0("node", nodes[node.mask])

  return( SRD.by.node)
}

plot.phylo(Trees$extant.sp.tree, no.margin = T, cex = 0.8)
nodelabels(cex = 0.7)

SRD.results.test <- SRD.Tree(tree = Trees$extant.sp.tree, 
                        mean.list = All.sp.data$means[mask.extant & mask.at.tree],
                        cov.matrix.list = All.sp.data$cov.mx[mask.extant & mask.at.tree],
                        sample.sizes = All.sp.data$n.sizes[mask.extant & mask.at.tree])

SRD.results.test$node135$SRD.result

#save(SRD.results, file = "Data/SRDatTreeResults.RData")

SRD.results$rolated <- SRD.results.test %>% ldply(., function(x) class(x$SRD.result) == "SRD") 

SRD.results$summary <- SRD.results.test[SRD.results$rolated[,2]] %>% laply(function(x) x$SRD.result$output[,5]) %>% alply(., 2, summary) %>% ldply(function (x) x) 
SRD.results$summary %>% melt %>%
  ggplot(., (aes(x=variable, y = value, group = X1))) +
  geom_label(aes(x=variable, y = value, group = X1, label = X1) )
  geom_point()


  SRD.node <- SRD.results.test[SRD.results$rolated[,2]]
> 


SRD.results.test [SRD.results$rolated[,2]] %>% llply(function(x) cbind(mean = mean(x$SRD.result$output[,1]), rs = x$RS.result[1]) )

SRD.selected <- list("Prosimian" = SRD.results.test$node71$SRD.result,
                     "Strepsirrhini" = SRD.results.test$node73$SRD.result,
                     "Lorisiformes" = SRD.results.test$node131$SRD.result,
                     "Galagidae" = SRD(ancestral.mx$Galago_senegalensis, ancestral.mx$Euoticus_elegantulus),
                     "Lorisidae" = SRD.results.test$node132$SRD.result,
                     "LemfxDaubs" = SRD.results.test$node74$SRD.result,
                     "Lemuriformes" = SRD.results.test$node75$SRD.result,
                     "LemxInd" = SRD.results.test$node99$SRD.result,
                     "Lemuridae" = SRD.results.test$node112$SRD.result,
                     "Indridae" = SRD.results.test$node100$SRD.result,
                     "LepxChe" = SRD.results.test$node76$SRD.result,
                     "Lepilemuridea" = SRD.results.test$node88$SRD.result,
                     "Cheirogaleidae" = SRD.results.test$node77$SRD.result
                     )


SRD.selected %>% llply(function (x) mean(x$output[,1]) )

SRD.select.genus <- list("Varecia" = SRD.results$node112,
                     "Lemur" = SRD.results$node125,
                     "Prolemur" = SRD.results$node126,
                     "Hapalemur" = SRD.results$node127,
                     "Eulemur" = SRD.results$node114,
                     "LemfxDaubs" = SRD.results$node74,
                     "Lemuriformes" = SRD.results$node75,
                     "LemxInd" = SRD.results$node99,
                     "Lemuridae" = SRD.results$node112,
                     "Indridae" = SRD.results$node100,
                     "LepxChe" = SRD.results$node76,
                     "Lepilemuridea" = SRD.results$node88,
                     "Cheirogaleidae" = SRD.results$node77
)


SRD.summary <- SRD.selected %>% laply(function(x) x$output[,5]) 
rownames(SRD.summary) <- names(SRD.selected) 


SRD.code <- SRD.selected %>% laply(function(x) x$model$code) 
rownames(SRD.code) <- names(SRD.selected) 
colnames(SRD.code) <- colnames(SRD.summary)

SRD.summary %<>% melt 
names(SRD.summary) <- c("node", "trait", "muDev")  
SRD.code %<>% melt
names(SRD.code) <- c("node", "trait", "code")  

SRD.to.plot <- cbind(SRD.summary, code =SRD.code$code)
names(SRD.to.plot)

ggplot(SRD.to.plot) +
  geom_point(aes(y = trait, x = muDev, size = code*0.4) , color = "darkgrey", alpha = 0.6) +
  geom_point(aes(y = trait, x = muDev, color = node), size = 2) +
    geom_vline(aes(xintercept=0))+
    facet_wrap(~node, ncol = 4) +  
  theme_bw() + panel_border(colour = "darkgrey", size = 01, remove = F) +
  background_grid(major = 'y', minor = "none") +
  theme(legend.position = "none", 
        strip.background = element_rect(fill = "transparent")) 


SRD.results.test$node77$SRD.result$model$code

SRD.sd.mean <- SRD.results.test[SRD.results$rolated[,2]] %>% ldply(., function(x) cbind(x$SRD.result$output, x$SRD.result$model$code, "A" = x$m.sd.cor[[1]], "B" = x$m.sd.cor[[2]][2:3]) ) 
names(SRD.sd.mean)[8] <- "code"
SRD.sd.mean $code


SRD.sd.mean$Oral <- rep(Modular.hyp[, "Oral"], 32)
SRD.sd.mean$Nasal <- rep(Modular.hyp[, "Nasal"], 32)
SRD.sd.mean$Zygomatic <- rep(Modular.hyp[, "Zygomatic"], 32)
SRD.sd.mean$Orbit <- rep(Modular.hyp[, "Orbit"], 32)
SRD.sd.mean$Nasal <- rep(Modular.hyp[, "Nasal"], 32)
SRD.sd.mean$Base <- rep(Modular.hyp[, "Base"], 32)
SRD.sd.mean$Vault <- rep(Modular.hyp[, "Vault"], 32)
SRD.sd.mean$Face <- rep(Modular.hyp[, "Face"], 32)
SRD.sd.mean$Neuro <- rep(Modular.hyp[, "Neuro"], 32)

SRD.sd.mean %>% dim()
1248/39


SRD.sd.mean %>% filter (.id == "node71"| 
                          .id == "node73"| 
                          .id == "node131"| 
                          .id == "node132"| 
                          .id == "node74"| 
                          .id == "node75"| 
                          .id == "node99"| 
                          .id == "node112"| 
                          .id == "node100"|
                          .id == "node76"| 
                          .id == "node88"| 
                          .id == "node77") %>%
  ggplot(.) +
  geom_point(aes(x = CSD, y = A.mean.cor, size = code, color = interaction(Oral, Nasal, Zygomatic, Orbit, Base, Vault) ), shape = 8) +
  #geom_smooth(aes(x = CSD, y = A.mean.cor, size = code), method = "lm", color = "red") +
  #geom_text(aes(x = CSD, y = A.mean.cor, label = A.trait, size = code), color = "red") +
  #geom_point(aes(x = CSD, y = A.sd.cor), color = "red", shape = 8) +
  geom_point(aes(x = CSD, y = B.mean.cor, size = code, color = interaction(Oral, Nasal, Zygomatic, Orbit, Base, Vault))) + 
  #geom_smooth(aes(x = CSD, y = B.mean.cor, size = code), color = "black", method = "lm") +
  #geom_text(aes(x = CSD, y = B.mean.cor, label = A.trait, size = code)) +
  #geom_point(aes(x = CSD, y = B.sd.cor), shape = 8) + 
  facet_wrap(~.id, ncol = 3)





SRD.sd.mean %>% filter (.id == "node71"| 
                          .id == "node73"| 
                          .id == "node131"| 
                          .id == "node132"| 
                          .id == "node74"| 
                          .id == "node75"| 
                          .id == "node99"| 
                          .id == "node112"| 
                          .id == "node100"|
                          .id == "node76"| 
                          .id == "node88"| 
                          .id == "node77") %>%
  ggplot(.) +
  geom_point(aes(x = CSD, y = A.sd.cor, size = code), color = "red") +
  geom_smooth(aes(x = CSD, y = A.sd.cor, size = code), method = "lm", color = "red") +
  #geom_text(aes(x = CSD, y = A.mean.cor, label = A.trait, size = code), color = "red") +
  #geom_point(aes(x = CSD, y = A.sd.cor), color = "red", shape = 8) +
  geom_point(aes(x = CSD, y = B.sd.cor, size = code)) + 
  geom_smooth(aes(x = CSD, y = B.sd.cor, size = code), method = "lm") +
  #geom_text(aes(x = CSD, y = B.mean.cor, label = A.trait, size = code)) +
  #geom_point(aes(x = CSD, y = B.sd.cor), shape = 8) + 
  facet_wrap(~.id, ncol = 3)

SRD.sd.mean %>% filter (.id == "node71"| 
                          .id == "node73"| 
                          .id == "node131"| 
                          .id == "node132"| 
                          .id == "node74"| 
                          .id == "node75"| 
                          .id == "node99"| 
                          .id == "node112"| 
                          .id == "node100"|
                          .id == "node76"| 
                          .id == "node88"| 
                          .id == "node77") %>%
  ggplot(.) +
  geom_point(aes(x = CSD, y = A.sd.cor/A.mean.cor), color = "red") +
  geom_text(aes(x = CSD, y = A.sd.cor/A.mean.cor, label = A.trait, size = code), color = "red") +
  geom_point(aes(x = CSD, y = B.sd.cor/B.mean.cor)) +
  geom_text(aes(x = CSD, y = B.sd.cor/B.mean.cor, label = A.trait, size = code)) +
  facet_wrap(~.id, ncol = 3)


Bvar.Wvar <- Drift.results$extant.sp$BW.compare %>% ldply(function(x) cbind(B.var = sum(diag(x$B.ed)), W.var = sum(diag(x$W)), r2W = CalcR2(x$W)) )
Bvar.Wvar %>% 
  filter (.id == "71"| 
            .id == "73"| 
            .id == "131"| 
            .id == "132"| 
            .id == "74"| 
            .id == "75"| 
            .id == "99"| 
            .id == "112"| 
            .id == "100"|
            .id == "76"| 
            .id == "88"| 
            .id == "77") %>%
  ggplot() +
  #geom_point(aes(x= r2W, y = B.var/W.var)) +
  geom_label(aes(x= r2W, y = B.var/W.var, label = .id))

plot.phylo(Trees$extant.sp.tree, no.margin = T, cex = 0.8)
nodelabels(text = round(Bvar.Wvar$B.var/Bvar.Wvar$W.var, digits=2), node = as.numeric(Bvar.Wvar$.id), adj = 0.8, bg = "grey")
nodelabels(text = round(Bvar.Wvar$r2W, digits=2), node = as.numeric(Bvar.Wvar$.id), adj = -0.3, bg = "yellow" )
legend("bottomleft", inset = 0.08,
       title ="",
       text.col = "grey10",
       legend = c("var.B/varW", "r2"), 
       fill = c("grey", "yellow") , 
       col = c("grey", "yellow") , 
       border = "grey", box.lwd = "n",
       bg= "transparent",
       cex = 2)


sum(diag(x$B.ed))/sum(diag(x$W))

CalcR2(Drift.results$extant.sp$BW.compare$`75`$W)
CalcR2(sp.main.data$Indri_indri$matrix$cov)



load("Data/Shapes.RData")
load("~/MScMadagascar/R/Shapes/file.RData")
load("~/MScMadagascar/R/Shapes/Aux.RData")
load("~/MScMadagascar/R/Shapes/shapes/traits.RData")
source('R/Shapes/ggshape.R', echo=TRUE)


SRD.plot.wire <- function (SRD.result, SHAPE, ROTACIONI =  c(1,-1,1), TTL = "Awesome SRD result") {
CODE <- SRD.result$model$code
muDev = as.numeric(SRD.result$output[,5])
as.numeric(muDev < -0.2)
as.numeric(muDev < 0)
COLORE <- abs(CODE * muDev)
COLORE <- CODE + as.numeric(muDev < -0.2) + as.numeric(muDev < 0)
SRD.relMean <- as.numeric(SRD.result$output[,1]) / mean(as.numeric(SRD.result$output[,1]))
plot.cheetows <- 
  ggshape(shape = SHAPE,
            wireframe = Aux $ tessel.39 [1:39, ],
            colors = rev(COLORE),
            rotation =ROTACIONI, 
            culo = 0.03, 
            thickness = 0.5) +
    #geom_point (aes (x = X, y = Y), alpha = 0.6, color ="darkgrey", size = 1.2) +
    #geom_label(aes (x = X, y = Y, label = traits ),  alpha = 0.4, size = 2, label.padding = unit(0.5, "mm")) +
    ggtitle(TTL) +
    theme(plot.title = element_text(face = "bold", size = 12),
          legend.position= "none", 
          legend.direction="horizontal",
          legend.text = element_text(size = 5),
          plot.margin = unit(c(0,0,-0,0), "cm"), 
          legend.key.size = unit(0.3, "cm"), 
          panel.margin.x = unit(0.3, "cm"), 
          panel.margin.y = unit(0.3, "cm") ) +
    scale_fill_gradientn(colors = c("#cccccc", "#fcae91", "#fb6a4a","#cb181d"),
                         guide = guide_colorbar(nbin=100, draw.ulim = T, draw.llim = T) ) +
    scale_color_gradientn(colors = c("#cccccc", "#fcae91", "#fb6a4a","#cb181d"),
                          #colors = myPalette(11),
                          guide = guide_colorbar(nbin=100, draw.ulim = T, draw.llim = T) )

avgSRD <- mean(as.numeric(SRD.result$output[,1]))

plot.doritows <- 
  ggshape(shape = SHAPE,
          wireframe = Aux $ tessel.39 [1:39, ],
          colors = rev(as.numeric(SRD.result$output[,1])),
          rotation =ROTACIONI, 
          culo = 0.03,
          thickness = 1) +
  #geom_point (aes (x = X, y = Y), alpha = 0.6, color ="darkgrey", size = 1.2) +
  #geom_label(aes (x = X, y = Y, label = traits ),  alpha = 0.4, size = 2, label.padding = unit(0.5, "mm")) +
  ggtitle(TTL) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position= "none", 
        legend.direction="horizontal",
        legend.text = element_text(size = 7),
        plot.margin = unit(c(0,0,-0,0), "cm"), 
        legend.key.size = unit(0.3, "cm"), 
        panel.margin.x = unit(0.3, "cm"), 
        panel.margin.y = unit(0.3, "cm") ) +
scale_fill_gradient2(high = "#92c5de", mid = "#cccccc", low ="#de2d26",  midpoint = avgSRD, labels = c(0, 0.2, 0.5, 0.8, 1), breaks = c(-0.1, 0.2, 0.5, 0.7, 1.1),
                    guide = guide_colorbar(nbin=1000, draw.ulim = T, draw.llim = T) ) +
scale_color_gradient2(high = "#92c5de", mid = "#cccccc", low = "#de2d26" ,  midpoint = avgSRD, labels = c(0, 0.2, 0.5, 0.8, 1), breaks = c(-0.1, 0.2, 0.5, 0.7, 1.1),
                     guide = guide_colorbar(nbin=10000, draw.ulim = T, draw.llim = T) )
#scale_fill_gradientn(colors = c("red", "darkgrey", "blue"),
#                       guide = guide_colorbar(nbin=100, draw.ulim = T, draw.llim = T) ) +
# scale_color_gradientn(colors = c("red", "darkgrey", "blue"),
#                       guide = guide_colorbar(nbin=100, draw.ulim = T, draw.llim = T) ) 
return(list(plot.muDev = plot.cheetows, 
            plot.relMean = plot.doritows , 
            color = COLORE, 
            SRD.relMean = SRD.relMean) ) 
}

SRD.plot.wire(SRD.result = SRD.results$node71, SHAPE = Shapes.sym$avahi, ROTACIONI =  c(-1,-1,1))

SRD.plot.wire(SRD.result = SRD.selected$Prosimian, SHAPE = Shapes.sym$tarsius, ROTACIONI =  c(-1,-1,1), TTL = "Prosimian\n Strepsirrhini + Tarsiidae")


SRD.selected.Tree.plot <- plot_grid(
  SRD.plot.wire(SRD.result = SRD.selected$Galagidae, SHAPE = Shapes.sym$galago, ROTACIONI =  c(-1,-1,1), TTL = "Galagidae")$plot.muDev,
  SRD.plot.wire(SRD.result = SRD.selected$Lorisidae, SHAPE = Shapes.sym$loris, ROTACIONI =  c(-1,-1,1), TTL = "Lorisidae")$plot.muDev,
  SRD.plot.wire(SRD.result = SRD.selected$Lorisiformes, SHAPE = Shapes.sym$nycticebus, ROTACIONI =  c(-1,-1,1), TTL = "Lorisiformes\n Galagidae x Lorisidae")$plot.muDev,
  
  SRD.plot.wire(SRD.result = SRD.selected$Lemuridae, SHAPE = Shapes.sym$varecia, ROTACIONI =  c(-1,-1,1), TTL = "Lemuridae")$plot.muDev,
  SRD.plot.wire(SRD.result = SRD.selected$Indridae, SHAPE = Shapes.sym$avahi, ROTACIONI =  c(-1,-1,1), TTL = "Indridae")$plot.muDev,
  SRD.plot.wire(SRD.result = SRD.selected$LemxInd, SHAPE = Shapes.sym$propithecus, ROTACIONI =  c(-1,-1,1), TTL = "Lemuridae x Indridae" )$plot.muDev,
  
  SRD.plot.wire(SRD.result = SRD.selected$Lepilemuridea, SHAPE = Shapes.sym$lepilemur, ROTACIONI =  c(1,-1,1), TTL = "Lepilemuridae")$plot.muDev,
  SRD.plot.wire(SRD.result = SRD.selected$Cheirogaleidae, SHAPE = Shapes.sym$microcebus, ROTACIONI =  c(1,-1,1), TTL = "Cheirogaleidae")$plot.muDev,
  SRD.plot.wire(SRD.result = SRD.selected$LepxChe, SHAPE = Shapes.sym$cheirogaleus, ROTACIONI =  c(1,-1,1), TTL = "Lepilemuridae x Cheirogaleidae")$plot.muDev,
  
  SRD.plot.wire(SRD.result = SRD.selected$Lemuriformes, SHAPE = Shapes.sym$hapalemur, ROTACIONI =  c(1,-1,1), TTL = "Lemuriformes")$plot.muDev,
  SRD.plot.wire(SRD.result = SRD.selected$LemfxDaubs, SHAPE = Shapes.sym$daubentonia, ROTACIONI =  c(1,-1,1), TTL = "Lemuriformes x Daubentonidae")$plot.muDev,
  SRD.plot.wire(SRD.result = SRD.selected$Strepsirrhini, SHAPE = Shapes.sym$perodicticus, ROTACIONI =  c(-1,-1,1), TTL = "Strepsirrhini")$plot.muDev, # )#,
  
  SRD.plot.wire(SRD.result = SRD.selected$Prosimian, SHAPE = Shapes.sym$tarsius, ROTACIONI =  c(-1,-1,1), TTL = "Prosimian\n Strepsirrhini + Tarsiidae" )$plot.muDev, 
  ncol = 3)

SRD.selected.Tree.plot


SRD.selected.Tree.plot <- plot_grid(
  SRD.plot.wire(SRD.result = SRD.selected$Galagidae, SHAPE = Shapes.sym$galago, ROTACIONI =  c(-1,-1,1), TTL = "Galagidae")$plot.relMean,
  SRD.plot.wire(SRD.result = SRD.selected$Lorisidae, SHAPE = Shapes.sym$loris, ROTACIONI =  c(-1,-1,1), TTL = "Lorisidae")$plot.relMean,
  SRD.plot.wire(SRD.result = SRD.selected$Lorisiformes, SHAPE = Shapes.sym$nycticebus, ROTACIONI =  c(-1,-1,1), TTL = "Lorisiformes\n Galagidae x Lorisidae")$plot.relMean,
  
  SRD.plot.wire(SRD.result = SRD.selected$Lemuridae, SHAPE = Shapes.sym$varecia, ROTACIONI =  c(-1,-1,1), TTL = "Lemuridae")$plot.relMean,
  SRD.plot.wire(SRD.result = SRD.selected$Indridae, SHAPE = Shapes.sym$avahi, ROTACIONI =  c(-1,-1,1), TTL = "Indridae")$plot.relMean,
  SRD.plot.wire(SRD.result = SRD.selected$LemxInd, SHAPE = Shapes.sym$propithecus, ROTACIONI =  c(-1,-1,1), TTL = "Lemuridae x Indridae" )$plot.relMean,
  
  SRD.plot.wire(SRD.result = SRD.selected$Lepilemuridea, SHAPE = Shapes.sym$lepilemur, ROTACIONI =  c(1,-1,1), TTL = "Lepilemuridae")$plot.relMean,
  SRD.plot.wire(SRD.result = SRD.selected$Cheirogaleidae, SHAPE = Shapes.sym$microcebus, ROTACIONI =  c(1,-1,1), TTL = "Cheirogaleidae")$plot.relMean,
  SRD.plot.wire(SRD.result = SRD.selected$LepxChe, SHAPE = Shapes.sym$cheirogaleus, ROTACIONI =  c(1,-1,1), TTL = "Lepilemuridae x Cheirogaleidae")$plot.relMean,
  
  SRD.plot.wire(SRD.result = SRD.selected$Lemuriformes, SHAPE = Shapes.sym$hapalemur, ROTACIONI =  c(1,-1,1), TTL = "Lemuriformes")$plot.relMean,
  SRD.plot.wire(SRD.result = SRD.selected$LemfxDaubs, SHAPE = Shapes.sym$daubentonia, ROTACIONI =  c(1,-1,1), TTL = "Lemuriformes x Daubentonidae")$plot.relMean,
  SRD.plot.wire(SRD.result = SRD.selected$Strepsirrhini, SHAPE = Shapes.sym$perodicticus, ROTACIONI =  c(-1,-1,1), TTL = "Strepsirrhini")$plot.relMean, 
  
  SRD.plot.wire(SRD.result = SRD.selected$Prosimian, SHAPE = Shapes.sym$tarsius, ROTACIONI =  c(-1,-1,1), TTL = "Prosimian\n Strepsirrhini + Tarsiidae" )$plot.relMean, 
  ncol = 3)

SRD.selected.Tree.plot


ggshape(shape = Shapes.sym$perodicticus,
        wireframe = Aux $ tessel.39 [1:39, ],
        colors = as.numeric(SRD.selected$Prosimian$output[,1]),
        rotation = c(-1,-1,1), 
        culo = 0.03,
        thickness = 1) +
  #geom_point (aes (x = X, y = Y), alpha = 0.6, color ="darkgrey", size = 1.2) +
  #geom_label(aes (x = X, y = Y, label = traits ),  alpha = 0.4, size = 2, label.padding = unit(0.5, "mm")) +
  ggtitle(TTL) +
  theme(plot.title = element_text(face = "bold", size = 12),
        #legend.position= "none", 
        legend.direction="horizontal",
        legend.text = element_text(size = 5),
        plot.margin = unit(c(0,0,-0,0), "cm"), 
        legend.key.size = unit(0.3, "cm"), 
        panel.margin.x = unit(0.3, "cm"), 
        panel.margin.y = unit(0.3, "cm") ) +
  #scale_fill_gradient(low = "#132B43", high = "#56B1F7", labels = c(0, 0.2, 0.5, 0.7, 1), breaks = c(0, 0.2, 0.5, 0.7, 1), 
  #                     guide = guide_colorbar(nbin=100, draw.ulim = T, draw.llim = T) ) +
  #scale_color_gradientn(colors = c("#7b3294", "#c2a5cf", "#cccccc", "#a6dba0","#008837"),
  #colors = myPalette(11),
  #                      guide = guide_colorbar(nbin=100, draw.ulim = T, draw.llim = T) )
  # scale_fill_gradientn(colors = c("red", "darkgrey", "blue"),
  #                       guide = guide_colorbar(nbin=100, draw.ulim = T, draw.llim = T) ) +
  # scale_color_gradientn(colors = c("red", "darkgrey", "blue"),
  #                       guide = guide_colorbar(nbin=100, draw.ulim = T, draw.llim = T) ) 

SRD.out <- SRD.selected %>% llply(function(x) x$output[, 6])

Mx.cor.mean <- Gen.cov.list %>% llply( function(x) data.frame( "trait" =  names(mean.sim(x)), "mean.cor" = mean.sim(abs(x)), "sd.cor" = mean.sd(abs(x)) ) , .progress = "text" )

Mx.cor.mean$W.Lorisiformes

SRD.selected$Prosimian$pc1
plot(SRD.selected$Prosimian)

SRD.selected.Tree.plot <- plot_grid(
  SRD.plot.wire(SRD.result = SRD(sp.main.data$Galago_senegalensis$matrix$cov, sp.main.data$Nycticebus_coucang$matrix$cov), SHAPE = Shapes.sym$galago, ROTACIONI =  c(-1,-1,1), TTL = "Galago x Nycticebus")$plot.muDev
  SRD.plot.wire(SRD.result = SRD.selected$Lorisidae, SHAPE = Shapes.sym$loris, ROTACIONI =  c(-1,-1,1), TTL = "Lorisidae")$plot.muDev,
  SRD.plot.wire(SRD.result = SRD.selected$Lorisiformes, SHAPE = Shapes.sym$nycticebus, ROTACIONI =  c(-1,-1,1), TTL = "Lorisiformes\n Galagidae x Lorisidae")$plot.muDev,
  SRD.plot.wire(SRD.result = SRD.selected$Lemuridae, SHAPE = Shapes.sym$varecia, ROTACIONI =  c(-1,-1,1), TTL = "Lemuridae")$plot.muDev,
  SRD.plot.wire(SRD.result = SRD.selected$Indridae, SHAPE = Shapes.sym$avahi, ROTACIONI =  c(-1,-1,1), TTL = "Indridae")$plot.muDev,
  SRD.plot.wire(SRD.result = SRD.selected$LemxInd, SHAPE = Shapes.sym$propithecus, ROTACIONI =  c(-1,-1,1), TTL = "Lemuridae x Indridae" )$plot.muDev,
  SRD.plot.wire(SRD.result = SRD.selected$Lepilemuridea, SHAPE = Shapes.sym$lepilemur, ROTACIONI =  c(1,-1,1), TTL = "Lepilemuridae")$plot.muDev,
  SRD.plot.wire(SRD.result = SRD.selected$Cheirogaleidae, SHAPE = Shapes.sym$microcebus, ROTACIONI =  c(1,-1,1), TTL = "Cheirogaleidae")$plot.muDev,
  SRD.plot.wire(SRD.result = SRD.selected$LepxChe, SHAPE = Shapes.sym$cheirogaleus, ROTACIONI =  c(1,-1,1), TTL = "Lepilemuridae x Cheirogaleidae")$plot.muDev,
  SRD.plot.wire(SRD.result = SRD.selected$Lemuriformes, SHAPE = Shapes.sym$hapalemur, ROTACIONI =  c(1,-1,1), TTL = "Lem-Ind x Lep-Che")$plot.muDev,
  SRD.plot.wire(SRD.result = SRD.selected$LemfxDaubs, SHAPE = Shapes.sym$daubentonia, ROTACIONI =  c(1,-1,1), TTL = "Lemuriformes\n Lemurs x Daubentonidae")$plot.muDev,
  SRD.plot.wire(SRD.result = SRD.selected$Strepsirrhini, SHAPE = Shapes.sym$perodicticus, ROTACIONI =  c(-1,-1,1), TTL = "Strepsirrhini")$plot.muDev, # )#,
  SRD.plot.wire(SRD.result = SRD.selected$Prosimian, SHAPE = Shapes.sym$tarsius, ROTACIONI =  c(-1,-1,1), TTL = "Prosimian\n Strepsirrhini + Tarsiidae" )$plot.muDev, 
  ncol = 3)

SRD.selected.Tree.plot


table(All.raw.main.data$All$Localidade == "ZOO")/2


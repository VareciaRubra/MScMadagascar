


SRD.Tree <- function (tree, cov.matrix.list, sample.sizes = NULL) 
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
                          tryCatch (expr = SRD(cov.matrices [[ indices [1] ]], cov.matrices [[ indices [2] ]]) , error = function(cond) vector("list") )
                        })
  names(SRD.by.node) <- paste0("node", nodes[node.mask])

  return( SRD.by.node)
}



plot.phylo(Trees$extant.sp.tree, no.margin = T, cex = 0.8)
nodelabels(cex = 0.7)

SRD.results <- SRD.Tree(tree = Trees$extant.sp.tree, 
                        cov.matrix.list = All.sp.data$cov.mx[mask.extant & mask.at.tree],
                        sample.sizes = All.sp.data$n.sizes[mask.extant & mask.at.tree])
#save(SRD.results, file = "Data/SRDatTreeResults.RData")

SRD.results$rolated <- SRD.results %>% ldply(., function(x) class(x) == "SRD") 

SRD.results$summary <- SRD.results[SRD.results$rolated[,2]] %>% laply(function(x) x$output[,5]) %>% alply(., 2, summary) %>% ldply(function (x) x) 

SRD.results %>% filter(names(SRD.results) == "node77" )
SRD.selected <- list("Prosimian" = SRD.results$node71,
                     "Strepsirrhini" = SRD.results$node73,
                     "Out.Madagascar" = SRD.results$node131,
                     "Madagascar" = SRD.results$node74,
                     "L-IxL-C" = SRD.results$node75,
                     "LemxInd" = SRD.results$node99,
                     "LepxChe" = SRD.results$node76)

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
    geom_point(aes(y = trait, x = muDev, color = node)) +
    geom_point(aes(y = trait, x = muDev, size = code*0.5) , color = "darkgrey", alpha = 0.6) +
    facet_wrap(~node) +  
  theme_bw() + panel_border(colour = "darkgrey", size = 01, remove = F) +
  background_grid(major = 'y', minor = "none") +
  theme(legend.position = "none", 
        strip.background = element_rect(fill = "transparent")) 


load("Data/Shapes.RData")
load("~/MScMadagascar/R/Shapes/file.RData")
load("~/MScMadagascar/R/Shapes/Aux.RData")
load("~/MScMadagascar/R/Shapes/shapes/traits.RData")
source('R/Shapes/ggshape.R', echo=TRUE)


SRD.plot.wire <- function (SRD.result, SHAPE, ROTACIONI =  c(1,-1,1), TTL = "Awesome SRD result") {
CODE <- SRD.result$model$code
muDev = as.numeric(SRD.result$output[,5])
plot.cheetows <- 
  ggshape(shape = SHAPE,
            wireframe = Aux $ tessel.39 [1:39, ],
            colors = CODE,
            rotation =ROTACIONI, 
            culo = 0.02, 
            thickness = 1) +
    geom_point (aes (x = X, y = Y), alpha = 0.6, color ="darkgrey", size = 1.2) +
    geom_label(aes (x = X, y = Y, label = traits ),  alpha = 0.6, size = 2, label.padding = unit(0.3, "mm")) +
    ggtitle(TTL) +
    theme(plot.title = element_text(face = "bold", size = 12),
          legend.position= "none", 
          legend.direction="horizontal",
          legend.text = element_text(size = 8),
          plot.margin = unit(c(0,0,-0,0), "cm"), 
          legend.key.size = unit(0.3, "cm"), 
          panel.margin.x = unit(0.3, "cm"), 
          panel.margin.y = unit(0.3, "cm") ) +
    scale_fill_gradientn(colors = c("darkgrey", "red"),
                         guide = guide_colorbar(nbin=100, draw.ulim = T, draw.llim = T) ) +
    scale_color_gradientn(colors = c("darkgrey", "red"),
                          #colors = myPalette(11),
                          guide = guide_colorbar(nbin=100, draw.ulim = T, draw.llim = T) )
  # scale_fill_gradientn(colors = c("red", "darkgrey", "blue"),
  #                       guide = guide_colorbar(nbin=100, draw.ulim = T, draw.llim = T) ) +
  # scale_color_gradientn(colors = c("red", "darkgrey", "blue"),
  #                       guide = guide_colorbar(nbin=100, draw.ulim = T, draw.llim = T) ) 
return(plot.cheetows)
}

SRD.selected.Tree.plot <- plot_grid(SRD.plot.wire(SRD.result = SRD.selected$LemxInd, SHAPE = Shapes.sym$propithecus, ROTACIONI =  c(-1,-1,1), TTL = "Lemuridae x Indridae" ),
SRD.plot.wire(SRD.result = SRD.selected$LepxChe, SHAPE = Shapes.sym$cheirogaleus, ROTACIONI =  c(1,-1,1), TTL = "Lepilemuridae x Cheirogaleidae"),
SRD.plot.wire(SRD.result = SRD.selected$`L-IxL-C`, SHAPE = Shapes.sym$hapalemur, ROTACIONI =  c(1,-1,1), TTL = "Lem-Ind x Lep-Che"),
SRD.plot.wire(SRD.result = SRD.selected$Madagascar, SHAPE = Shapes.sym$daubentonia, ROTACIONI =  c(1,-1,1), TTL = "Madagascar\n Lemurs x Daubentonidae"),
SRD.plot.wire(SRD.result = SRD.selected$Out.Madagascar, SHAPE = Shapes.sym$loris, ROTACIONI =  c(-1,-1,1), TTL = "Out Madagascar\n Galagidae x Lorisidae"),
SRD.plot.wire(SRD.result = SRD.selected$Strepsirrhini, SHAPE = Shapes.sym$galago, ROTACIONI =  c(-1,-1,1), TTL = "Strepsirrhini"),
SRD.plot.wire(SRD.result = SRD.selected$Prosimian, SHAPE = Shapes.sym$tupaia, ROTACIONI =  c(-1,-1,1), TTL = "Prosimian\n Strepsirrhini + Tarsiidae" ))

SRD.selected.Tree.plot

as.numeric(SRD.selected$Prosimian$output[,5])

plot.phylo(Trees$extant.sp.tree, no.margin = T, cex = 0.8)
nodelabels(cex = 0.7)
plot(SRD.results$SRD.by.node$'100')
plot(SRD(mx.list.taxonomy$W.Indridae, mx.list.taxonomy$W.Lemuridae))

plot(SRD.results$SRD.by.node$'76')
plot(SRD(mx.list.taxonomy$W.Cheirogaleidae, mx.list.taxonomy$Lepilemur))

plot(SRD(cov.matrices[as.character(tree$edge[,2][tree$edge[,1] == 100])])[1,2][[1]] )


cov.matrices[as.character(tree$edge[,2][tree$edge[,1] == as.character(100)]) ]



as.character(tree$edge[,2][tree$edge[,1] == as.character(100)]) 

head(SRD.code)

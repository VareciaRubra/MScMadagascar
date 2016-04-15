


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
SRD.results$rolated <- SRD.results %>% ldply(., function(x) class(x) == "SRD") 

SRD.results$summary <- SRD.results[SRD.results$rolated[,2]] %>% laply(function(x) x$output[,5]) %>% alply(., 2, summary) %>% ldply(function (x) x) 

SRD.results %>% filter(names(SRD.results) == "node77" )
SRD.selected <- list("Prosimian" = SRD.results$node71,
                     "Strepsirrhini" = SRD.results$node73,
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


Shaper <- function( SHAPE){
  SHAPE.lm <- SHAPE [, 1]
  SHAPE.skull <- array (as.matrix (SHAPE [, 2:7]), c(27, 2, 3, 2))
  SHAPE.skull <- aperm(SHAPE.skull, c(1, 3, 2, 4))
  SHAPE.left <- SHAPE.skull [, , 1, ]
  SHAPE.right <- SHAPE.skull [, , 2, ]
  
  SHAPE.lm <- as.character (SHAPE.lm)
  
  SHAPE.lm[is.na(SHAPE.lm)] <-  'NA'
  SHAPE.lm %<>% gsub("La", "LA", .)
  SHAPE.lm %<>% gsub("Zi", "ZI", .)
  SHAPE.lm %<>% gsub("As", "AS", .)
  SHAPE.lm %<>% gsub("Ts", "TS", .)
  SHAPE.lm %<>% gsub("Oc", "OC", .)
  
  dimnames (SHAPE.left) <- list (SHAPE.lm[1:27], c ('X', 'Y', 'Z'))
  dimnames (SHAPE.right) <- list (SHAPE.lm[28:54], c ('X', 'Y', 'Z'))
  
  
  SHAPE.full <- glue.skulls(SHAPE.left, SHAPE.right) [[1]]
  
  SHAPE.av <- procGPA(SHAPE.full) $ mshape
  
  dimnames (SHAPE.av) <- dimnames (SHAPE.full) [1:2]
  rownames (SHAPE.full)
  
  SHAPE.sym <- OSymm (SHAPE.av, 1:9, 28:45, 10:27) $ symmconf
  
  rownames (SHAPE.sym) <- gsub ('d', '-D', gsub ('e', '-E', rownames (SHAPE.sym)))
  
  
  SHAPE.sym <- SHAPE.sym [which (rownames (SHAPE.sym) %in% rownames (hapa.sym)), ]
  
  SHAPE.sym <- SHAPE.sym [match(rownames (hapa.sym), rownames (SHAPE.sym)), ]
  
  return(SHAPE.sym)
}



CODE <- SRD.results$node71$model$code
muDEV <- SRD.results$node71$output[,5]
  myBreaks <- c(round(min(muDEV), digits = 1), round(mean (muDEV), digits = 1 ), round(max(muDEV), digits = 1 )) 
 
    ggshape(shape = ,
            wireframe = Aux $ tessel.39 [1:39, ],
            colors = muDEV,
            rotation = ROTACIONI , 
            culo = 0.02, 
            thickness = 1) +
    geom_point (aes (x = X, y = Y), alpha = 0.6, color ="darkgrey", size = 1.2) +
    geom_label(aes (x = X, y = Y, label = traits ),  alpha = 0.6, size = 3, label.padding = unit(0.3, "mm")) +
    #ggtitle(paste("PC", i, sep = " ")) +
    theme(plot.title = element_text(face = "bold", size = 12),
          legend.position= c(0.3,0.1), 
          legend.direction="horizontal",
          legend.text = element_text(size = 8),
          plot.margin = unit(c(0,0,-0,0), "cm"), 
          legend.key.size = unit(0.3, "cm"), 
          panel.margin.x = unit(0.3, "cm"), 
          panel.margin.y = unit(0.3, "cm") ) +
    scale_fill_gradientn(paste("PC", i, sep = " "), limits = c(min(WPCs[, i])-0.09, max(WPCs[, i])+0.09),
                         breaks= myBreaks,
                         #colors = myPalette(11),
                         colors = myPalette(11),
                         guide = guide_colorbar(nbin=100, draw.ulim = T, draw.llim = T) ) +
    scale_color_gradientn(paste("PC", i, sep = " "), limits = c(min(WPCs[, i]) -0.09, max(WPCs[, i]) +0.09 ),
                          breaks= myBreaks,
                          #colors = myPalette(11),
                          colors = myPalette(11),
                          guide = guide_colorbar(nbin=100, draw.ulim = T, draw.llim = T) ) 
  


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

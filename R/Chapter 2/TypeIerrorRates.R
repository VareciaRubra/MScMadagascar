DriftTestInterpretation <- function (tree, mean.list, cov.matrix.list, W.fixed = NULL, sample.sizes = NULL) 
{
  if (!all(tree$tip.label %in% names(mean.list))) 
    stop("All tip labels must be in names(mean.list).")
  if (!all(tree$tip.label %in% names(cov.matrix.list))) 
    stop("All tip labels must be in names(cov.matrix.list).")
  cov.matrices <- PhyloW(tree, cov.matrix.list, sample.sizes)
  
    if(is.null(W.fixed)) W.fixed <- cov.matrices[[length(cov.matrices)]] else W.fixed <-W.fixed
  
  nodes <- names(cov.matrices)
  node.mask <- laply(nodes, function(x) length(getMeans(mean.list, tree, x))) > 3
  if (!any(node.mask)) 
    stop("For the regression test:/nAt least one node must have more than 4 descendents in mean.list")
  
  ind.cont<- ldply(mean.list, function(x) x) 
  rownames(ind.cont) <- ind.cont[,1]
  ind.cont <- ind.cont[,-1]
  ind.cont <- apply(ind.cont, 2, FUN = function (x) ape::pic(x, tree) ) 
  
  BW.compare <- llply(nodes[node.mask], function(node) DumBW.compare(means = getMeans(mean.list, tree, node), 
                                                                     contrasts = getContrasts(ind.cont, tree, node), 
                                                                     W.mx = cov.matrices[[node]])
  )
  names(BW.compare) <- nodes[node.mask]
  test.list.cor <- llply(nodes[node.mask], function(node) PCScoreCorr(means = getMeans(mean.list, tree, node), 
                                                                      taxons = names(getMeans(mean.list, tree, node)), 
                                                                      cov.matrix =  W.fixed, 
                                                                      show.plots = FALSE,
                                                                      title.plot =  node) )
  names(test.list.cor) <- nodes[node.mask]
  
  test.list.reg <- llply(nodes[node.mask], function(node) DriftTest0(means = getMeans(mean.list, tree, node), 
                                                                               cov.matrix = W.fixed, 
                                                                               show.plot = FALSE))
  names(test.list.reg) <- nodes[node.mask]
  
  return(list ("Correlation.W.fixed" = test.list.cor,
               "Regression.W.fixed" = test.list.reg,
               "BW.compare" = BW.compare
) )
}

Drift.results$fixed <- DriftTestInterpretation(tree = Trees$all.with.ed, 
                        mean.list = All.sp.data$means[mask.at.tree][rowSums(missing.ed) == 39], 
                        cov.matrix.list = All.sp.data$cov.mx[mask.at.tree][rowSums(missing.ed) == 39], 
                        W.fixed = NULL,
                        #W.fixed = ancestral.mx$`42`,
                        sample.sizes = All.sp.data$n.sizes[mask.at.tree][rowSums(missing.ed) == 39] )

Drift.results$fixed <- DriftTestInterpretation(tree = Trees$extant.sp.tree, 
                                               mean.list = All.sp.data$means[mask.extant & mask.at.tree], 
                                               cov.matrix.list = All.sp.data$cov.mx[mask.extant & mask.at.tree], 
                                               W.fixed = NULL,
                                               sample.sizes = All.sp.data$n.sizes[mask.extant & mask.at.tree] )
Drift.results$fixed <- DriftTestInterpretation(tree = Trees$extant.sp.tree, 
                                               mean.list = All.sp.data$means[mask.extant & mask.at.tree], 
                                               cov.matrix.list = All.sp.data$cov.mx[mask.extant & mask.at.tree], 
                                               W.fixed = Gen.cov.list$Lepilemur,
                                               sample.sizes = All.sp.data$n.sizes[mask.extant & mask.at.tree] )

Drift.results.Toplot$Fixed$Plots$Corr.Contrasts  <- Drift.results$fixed$Correlation.W.fixed %>% llply(function (x) x$P.value.plot) %>% llply(function(x) x + theme(legend.position = "none") ) %>% cowplot::plot_grid(plotlist = .) 
Drift.results.Toplot$FixedW.Extants$Results$Node.ref <- Drift.results$fixed$Correlation.W.fixed %>% ldply(function(x) dim(table(x$Bonferroni == T, exclude = T)) ==2) %>% .[,1]
Drift.results.Toplot$FixedW.Extants$Results$Corr.Ed.1 <- Drift.results$fixed$Correlation.W.fixed  %>% ldply(function(x) (table(x$Bonferroni == T) %>% dimnames() %>% .[[1]] %>% length()) ==2)  %>% .[,2]
Drift.results.Toplot$FixedW.Extants$Results$Regression.Ed <- Drift.results$fixed$Regression.W.fixed %>% ldply(function(x) x$drift_rejected )  %>% .[,2]

Plot.Drift.Results(tree = Trees$all.with.ed, 
                   info = Drift.results.Toplot$All.sp$Results, 
                   font.size = 0.7,
                   contrasts = F,
                   paleta.deriva = c("#ECA400", "#1BAAB7", "#CBC9E0"))

Plot.Drift.Results(tree = Trees$extant.sp.tree, 
                   info = Drift.results.Toplot$FixedW.Extants$Results, 
                   font.size = 0.6,
                   contrasts = F,
                   paleta.deriva = c("#ECA400", "#1BAAB7", "#CBC9E0"))



Drift.results.Toplot$Fixed$Plots$Corr.Ed.1  <- Drift.results$fixed$Correlation.W.fixed %>% llply(function (x) x$P.value.plot) %>% llply(function(x) x + theme(legend.position = "none") ) %>% cowplot::plot_grid(plotlist = .) 

Drift.results$fixed$Correlation.W.fixed$`135`$P.value.plot

Correlation.Tree.plot.Wfixed <- plot_grid(
  Drift.results$fixed$Correlation.W.fixed$`131`$P.value.plot + theme(legend.position = "none") + ggtitle ("Lorisiformes") ,
  Drift.results$fixed$Correlation.W.fixed$`132`$P.value.plot + theme(legend.position = "none") + ggtitle ("Lorisidae") ,
  Drift.results$fixed$Correlation.W.fixed$`134`$P.value.plot + theme(legend.position = "none") + ggtitle ("Nyc-Lor"), 
  
  Drift.results$fixed$Correlation.W.fixed$`99`$P.value.plot + theme(legend.position = "none") + ggtitle ("Lemuridae x Indriidae") ,
  Drift.results$fixed$Correlation.W.fixed$`100`$P.value.plot + theme(legend.position = "none") + ggtitle ("Indriidae") ,
  Drift.results$fixed$Correlation.W.fixed$`112`$P.value.plot + theme(legend.position = "none") + ggtitle ("Lemuridae") ,
  
  
  Drift.results$fixed$Correlation.W.fixed$`76`$P.value.plot + theme(legend.position = "none") + ggtitle ("Lepilemuridae x Cheirogaleidae") ,
  Drift.results$fixed$Correlation.W.fixed$`77`$P.value.plot + theme(legend.position = "none") + ggtitle ("Cheirogaleidae") ,
  Drift.results$fixed$Correlation.W.fixed$`88`$P.value.plot + theme(legend.position = "none") + ggtitle ("Lepilemuridae") ,
  
  Drift.results$fixed$Correlation.W.fixed$`73`$P.value.plot + theme(legend.position = "none") + ggtitle ("Strepsirrhini") ,
  Drift.results$fixed$Correlation.W.fixed$`74`$P.value.plot + theme(legend.position = "none") + ggtitle ("Lemuriformes x Daubentonidae") ,
  Drift.results$fixed$Correlation.W.fixed$`75`$P.value.plot + theme(legend.position = "none") + ggtitle ("Lemuriformes") ,
  Drift.results$fixed$Correlation.W.fixed$`71`$P.value.plot + theme(legend.position = "none") + ggtitle ("Prosimian") ,
  ncol = 3)

Correlation.Tree.plot.Wfixed
save_plot(filename = "Figures/Correlation_Tree_plot_Wfixed.pdf", plot = Correlation.Tree.plot.Wfixed, 
          base_aspect_ratio = 0.3, base_height = 15, base_width = 9.5)


Regression.Tree.plot.fixed.W <- Drift.results$extant.sp$Regression.test %>% llply(Plot.Drift.regression)

Regression.Tree.plot.W.fixed <- plot_grid(
  Regression.Tree.plot.fixed.W$`134`$plot  + ggtitle ("Nyc-Lor"), 
  Regression.Tree.plot.fixed.W$`132`$plot  + ggtitle ("Lorisidae") ,
  Regression.Tree.plot.fixed.W$`131`$plot + ggtitle ("Lorisiformes") ,
  
  Regression.Tree.plot.fixed.W$`112`$plot + ggtitle ("Lemuridae") ,
  Regression.Tree.plot.fixed.W$`100`$plot  + ggtitle ("Indriidae") ,
  Regression.Tree.plot.fixed.W$`99`$plot  + ggtitle ("Lemuridae x Indriidae") ,
  
  Regression.Tree.plot.fixed.W$`88`$plot  + ggtitle ("Lepilemuridae") ,
  Regression.Tree.plot.fixed.W$`77`$plot + ggtitle ("Cheirogaleidae") ,
  Regression.Tree.plot.fixed.W$`76`$plot + ggtitle ("Lepilemuridae x Cheirogaleidae") ,
  
  Regression.Tree.plot.fixed.W$`75`$plot  + ggtitle ("Lemuriformes") ,
  Regression.Tree.plot.fixed.W$`74`$plot  + ggtitle ("Lemuriformes x Daubentonidae") ,
  Regression.Tree.plot.fixed.W$`73`$plot  + ggtitle ("Strepsirrhini") ,
  
  Regression.Tree.plot.fixed.W$`71`$plot  + ggtitle ("Prosimian") ,
  ncol = 3)

Regression.Tree.plot.W.fixed
save_plot(filename = "Figures/Regression_Tree_plot_Wfixed.pdf", plot = Regression.Tree.plot.W.fixed, 
          base_aspect_ratio = 0.3, base_height = 15, base_width = 9.5)

tabelinha <- Regression.Tree.plot.fixed.W  %>% ldply(function (x) x$tabelinha)
rownames(tabelinha) <- tabelinha[,1]

tabelinha[-1] %>% xtable()








GetBW<- function (tree, mean.list, cov.matrix.list, sample.sizes = NULL) 
{
  if (!all(tree$tip.label %in% names(mean.list))) 
    stop("All tip labels must be in names(mean.list).")
  if (!all(tree$tip.label %in% names(cov.matrix.list))) 
    stop("All tip labels must be in names(cov.matrix.list).")
  cov.matrices <- PhyloW(tree, cov.matrix.list, sample.sizes)
  
  nodes <- names(cov.matrices)
  node.mask <- laply(nodes, function(x) length(getMeans(mean.list, tree, x))) > 3
  if (!any(node.mask)) 
    stop("For the regression test:/nAt least one node must have more than 4 descendents in mean.list")
  
  ind.cont<- ldply(mean.list, function(x) x) 
  rownames(ind.cont) <- ind.cont[,1]
  ind.cont <- ind.cont[,-1]
  ind.cont <- apply(ind.cont, 2, FUN = function (x) ape::pic(x, tree) ) 
  
  BW.compare <- llply(nodes[node.mask], function(node) DumBW.compare(means = getMeans(mean.list, tree, node), 
                                                                     contrasts = getContrasts(ind.cont, tree, node), 
                                                                     W.mx = cov.matrices[[node]])  )
  names(BW.compare) <- nodes[node.mask]
  
  test.list.reg <- llply(nodes[node.mask], function(node) DriftTest0(means = getMeans(mean.list, tree, node), 
                                                                     cov.matrix = cov.matrices[[node]], 
                                                                     show.plot = FALSE))
  names(test.list.reg) <- nodes[node.mask]
  
  return(list ("Observed.Regression.Test" = test.list.reg,
               "BW.compare" = BW.compare
  ) )
}

BWtoSimulate <- vector("list")

BWtoSimulate$All.sp <- GetBW(tree = Trees$all.with.ed, 
                             mean.list = All.sp.data$means[mask.at.tree][rowSums(missing.ed) == 39], 
                             cov.matrix.list = All.sp.data$cov.mx[mask.at.tree][rowSums(missing.ed) == 39], 
                             sample.sizes = All.sp.data$n.sizes[mask.at.tree][rowSums(missing.ed) == 39] )

BWtoSimulate$Extants <- GetBW(tree = Trees$extant.sp.tree, 
                              mean.list = All.sp.data$means[mask.extant & mask.at.tree], 
                              cov.matrix.list = All.sp.data$cov.mx[mask.extant & mask.at.tree], 
                              sample.sizes = All.sp.data$n.sizes[mask.extant & mask.at.tree] )

var.B <- apply(BWtoSimulate$Extants$BW.compare$`76`$Proj.B, 2, var)
var.B/sum(var.B) 



BWtoSimulate$Extants$BW.compare %>% llply( function (x) apply(x$Proj.B, 2, var) )

W.percent <- BWtoSimulate$Extants$BW.compare%>% ldply(function (x) eigen(x$W)$values[1:10]/sum(eigen(x$W)$values))
names(W.percent) <- c(".node", paste0("PC", 1:10))
B.percent <- BWtoSimulate$Extants$BW.compare %>% ldply( function (x) apply(x$Proj.B[,1:10], 2, var)/ sum(apply(x$Proj.B[,1:10], 2, var)) )
names(B.percent) <- c(".node", paste0("PC", 1:10))
W.fixedpcpercent <- eigen(Gen.cov.list$W.Prosimian)$values[1:10]/sum(eigen(Gen.cov.list$W.Prosimian)$values)

B.percent <- BWtoSimulate$Extants$BW.compare %>% ldply( function (x) apply(x$Proj.B, 2, var) )
names(B.percent) <- c(".node", paste0("PC", 1:39))
W.fixedpcpercent <- eigen(Gen.cov.list$W.Prosimian)$values


#W.fixedpcpercent <- as.data.frame(W.fixedpcpercent)
#rownames(W.fixedpcpercent)<- paste0("PC", 1:10)
#W.fixedpcpercent$PC <- rownames(W.fixedpcpercent)
#names(W.fixedpcpercent) <- c("value", "PC")

B.percent$Matrix <- rep("B", dim(B.percent)[1])
W.percent$Matrix <- rep("W", dim(W.percent)[1])

B.percent %>% names
W.percent %>% names
W.fixedpcpercent.to.be <- B.percent
for (i in 1:  dim(W.fixedpcpercent.to.be)[1]) W.fixedpcpercent.to.be[i,2:40] <- W.fixedpcpercent
W.fixedpcpercent.to.be$Matrix <- "W.fixed"

pcpercent.pernode <- rbind(B.percent,# W.percent, 
                           W.fixedpcpercent.to.be)



pcpercent.pernode.selected <- pcpercent.pernode %>% melt %>% 
  filter (.node == "71"| 
            .node == "75"|
            .node == "74"| 
            .node == "73"| 
            .node == "88"| 
            .node == "77"|
            .node == "76"| 
            .node == "100"| 
            .node == "112"| 
            .node == "99"|
            .node == "134"| 
            .node == "132"| 
            .node == "131" ) 

pcpercent.pernode.selected$.node <- factor(pcpercent.pernode.selected$.node, levels = rev(c("71", 
                                                                                            "75", "74", "73", 
                                                                                            "88", "77", "76",
                                                                                            "112", "100", "99", 
                                                                                            "134", "132" , "131" ) ) )
                                                                                        
pcpercent.pernode.selected$variable %>% class                                                                                       

pcpercent.pernode.selected %>%
  ggplot( aes( x = variable, y = log(value), group = Matrix)  )+
    geom_point(aes(shape = Matrix, color = Matrix), alpha = 0.5, size = 0.5) +
    geom_line(aes(group = Matrix, color = Matrix), alpha = 0.4) +
  scale_color_manual(values = c("black",  "red") )+
    facet_wrap(~.node, ncol = 3) +
    theme_bw() + theme_minimal() + theme(axis.text.x = element_text(angle = 270, size = 7)) +
 ylab("log variance in each PC") + xlab("")



P <- BWtoSimulate$Extants$BW.compare$`71`$B.ed
G <- BWtoSimulate$Extants$BW.compare$`71`$W
pop <- BWtoSimulate$Extants$BW.compare$`71`$B.sample.size



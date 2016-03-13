require (geomorph)
require (shapes)
require (RColorBrewer)
require (ape)
require (Morphometrics)
require (expm)
require (plyr)
require (plotrix)
require (doMC)
registerDoMC (cores = 6)
require (reshape2)
require (ggplot2)
require (phytools)
require (geiger)
require (mvtnorm)
require (psych)
require (MCMCglmm)
require (boa)
require (grid)
require (gridExtra)
require (gridBase)
require (surface)
require (rmarkdown)
require (knitr)
require (pander)
require (shape)
require (OUwie)

## require (devtools)
## install_github('uyedaj/bayou')
## install_github ('mkoohafkan/kfigr')

attach ('../../Databases/Reference.RData')
attach ('../../Databases/ED.RData')
attach ('../../Databases/OneDef.RData')
attach ('../../Databases/Tree.RData')
attach ('../../Databases/Aux.RData')
attach ('../../Databases/LifeHistory.RData')
#attach ('../../covTensor/Work/post.vcv.RData')

options(contrasts = c('contr.sum', 'contr.poly'))

.source.files <- dir('../FuncR/', pattern = '.R', full.names = TRUE)
.source.files <- .source.files [!grepl ('~', .source.files)]
for (i in 1:length (.source.files))
  source (.source.files [i])

## Tree $ NodeNames [1, 2] <- 'Anthropoidea'
## save(Tree, file = '../../Databases/Tree.RData')
## detach (file:../../Databases/Tree.RData)

## OneDef <-
##   llply (OneDef, function (L)
##          {
##            colnames (L $ local) <- rownames (Aux $ def.hyp)
##            dimnames (L $ ml.vcv) <- list (rownames (Aux $ def.hyp),
##                                           rownames (Aux $ def.hyp))
##            colnames (L $ group.mean) <- rownames (Aux $ def.hyp)
##            names (L $ mean) <- rownames (Aux $ def.hyp)
##            L
##          })

## save (OneDef, file = '../../Databases/OneDef.RData')
## detach (file:../../Databases/OneDef.RData)
## rm (OneDef)

allo.Data <- list ()
allo.Plots <- list ()

allo.Data <-
  within (allo.Data,
          {

            BW.Def <- 
              alply (110:217, 1,
                     MLMem, tree = Tree [[1]],
                     data = OneDef, what = 'local',
                     .parallel = TRUE)

            BW.ED <- 
              alply (110:217, 1,
                     MLMem, tree = Tree [[1]],
                     df = ldply (ED, function (L) L $ ed[, -20]),
                     .parallel = TRUE)

            names (BW.Def) <- Tree $ NodeNames [, 2]
            names (BW.ED) <- names (BW.Def)

            W.node.ED <- llply (BW.ED, function (L) L $ W)
            names (W.node.ED) <- names (BW.ED)
            
            W.all.ED <- c(llply (ED, function (L) L $ ed.vcv [-20, -20]), W.node.ED)
            names (W.all.ED) <- c(names (ED), names (BW.Def))
            
            B.ED <- llply(BW.ED, function (L) L $ B)
            names (B.ED) <- names (BW.ED)

            mean.ED <- llply(BW.ED, function (L) L $ means)

            W.node.Def <- llply (BW.Def, function (L) L $ W)
            names (W.node.Def) <- names (BW.Def)
            
            W.all.Def <- c(llply (OneDef, function (L) L $ ml.vcv), W.node.Def)
            names (W.all.Def) <- c(names (OneDef), names (BW.Def))
            
            B.Def <- llply(BW.Def, function (L) L $ B)
            names (B.Def) <- names (BW.Def)

            mean.Def <- llply(BW.Def, function (L) L $ means)

            df.res <- laply(BW.Def, function (L) L $ df.res)
          })

## CAC RSC Ancestral
allo.Data <-
  within (allo.Data,
          {
            n.node <- 108

            n.tip <- laply(mean.Def, nrow)

            wAC.node.Def <- laply (W.node.Def, function (C) c(C [-1, 1] / C [1, 1]))
            rownames(wAC.node.Def) <- Tree $ NodeNames [, 2]

            ### Directions
            root.CAC <- Normalize (wAC.node.Def [1, ])
            names (root.CAC) <- rownames (Aux $ def.hyp) [-1]
            
            root.Iaa <- diag (38) - (root.CAC %*% t(root.CAC))

            root.RSC <-
              eigen (t (root.Iaa) %*% W.node.Def [[1]] [-1, -1] %*% root.Iaa) $ vectors
            colnames (root.RSC) <- paste ('RSC', 1:38, sep = '')
            rownames (root.RSC) <- rownames (W.node.Def [[1]]) [-1]

            ### Arrangement
            families <- Tree $ NodeNames [grep ('idae', Tree $ NodeNames [, 2]), ]
            families <- families [!grepl ('\\+', families [, 2]), ]
            subfamilies <- Tree $ NodeNames [grep ('inae', Tree $ NodeNames [, 2]), ]
            tribes <- Tree $ NodeNames [grep ('ini', Tree $ NodeNames [, 2]), ] [-(1:3), ]
            taxo.arr <- rbind (families [-8,], subfamilies [3,], tribes)

            rm (families, subfamilies, tribes)
            
            raw.df <- ldply (OneDef, function (L) L $ local)

            taxo.alleg <- c()

            for (i in 1:nrow (taxo.arr))
              {
                subtree <- extract.clade (Tree [[1]], taxo.arr [i, 1])
                taxo.alleg [raw.df $ .id %in% subtree $ tip.label] <- taxo.arr [i, 2]
              }

            rm (subtree)
          
            onedef.df <-
              cbind (raw.df [, 1:2],
                     as.matrix (raw.df [, -(1:2)]) %*% as.matrix (root.CAC),
                     as.matrix (raw.df [, -(1:2)]) %*% root.RSC [, 1:2])

            rm (raw.df, taxo.arr)

            onedef.df $ group <- taxo.alleg
            onedef.df <- onedef.df [, c(1, 6, 2, 3, 4, 5)]
            
            colnames (onedef.df) <- c('animal', 'taxo.group', 'logCS',
                                      'CACw', 'RSC1w', 'RSC2w')
            logCS.mean <- mean (onedef.df $ logCS)
            ### CENTER
            onedef.df $ logCS <- onedef.df $ logCS - logCS.mean
            
            onedef.mean <- ddply (onedef.df, .(animal), summarize,
                                  mCS = mean (logCS),
                                  mCAC = mean (CACw),
                                  mRSC1 = mean (RSC1w),
                                  mRSC2 = mean (RSC2w))
            
            rm (taxo.alleg)
          })

### Model
set.seed (42)

allo.Data <-
  within (allo.Data,
          {
            
            allo.mcmc.randomreg <-
              MCMCglmm(CACw ~ logCS,
                       random = ~ us(1 + logCS):animal,
                       rcov = ~ units,
                       data = onedef.df, pr = TRUE,
                       nitt = 100000, burnin = 50000, thin = 50,
                       pedigree = Tree [[1]],
                       family = 'gaussian')
          })

summary (allo.Data $ allo.mcmc.randomreg) 

#plot (allo.Data $ allo.mcmc.randomreg)

allo.Data <-
  within (allo.Data,
          {
            allo.mcmc.randomreg.idh <-
              MCMCglmm(CACw ~ logCS,
                       random = ~ idh(1 + logCS):animal,
                       rcov = ~ units,
                       data = onedef.df, pr = TRUE,
                       nitt = 100000, burnin = 50000, thin = 50,
                       pedigree = Tree [[1]],
                       family = 'gaussian')
          })

summary (allo.Data $ allo.mcmc.randomreg.idh) 

#plot (allo.Data $ allo.mcmc.randomreg.idh)


### Plot Preparation 1: Scatterplot
allo.Data <-
  within (allo.Data,
          {
            allo.mcmc.randomreg.us <- allo.mcmc.randomreg

            allo.mcmc.randomreg <- allo.mcmc.randomreg.idh
            
            allo.CAC.hull <- ddply (onedef.df, .(animal), summarize,
                                    hCS = logCS [chull(logCS, CACw)],
                                    hCAC = CACw [chull(logCS, CACw)])
            
            allo.CAC.group.hull <- ddply (onedef.df, .(taxo.group), summarize,
                                          hCS = logCS [chull(logCS, CACw)],
                                          hCAC = CACw [chull(logCS, CACw)])
            
          })

allo.Plots $ CAC.vs.logCS <- 
ggplot(allo.Data $ onedef.mean) +
  geom_point(aes(x = logCS, y = CACw, colour = taxo.group),
             allo.Data $ onedef.df, shape = '+') +
  geom_polygon(aes (x = hCS, y = hCAC, group = animal),
               allo.Data $ allo.CAC.hull, alpha = 0.05) +
  geom_polygon(mapping = aes (x = hCS, y = hCAC, colour = taxo.group,
                 fill = taxo.group),
               allo.Data $ allo.CAC.group.hull, alpha = 0.3) +
  theme_bw() +
  xlab ('log Centroid Size') +
  ylab ('Common Allometric Component') +
  scale_color_brewer('Clade', palette = 'Paired') +
  scale_fill_brewer('Clade', palette = 'Paired') +
  theme(legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit (0.25, 'cm'),
        legend.position = c(0.08,0.8),
        legend.background = element_rect(colour='black', size = 0.05))

### Plot Preparation 2: Phylogeny W Slopes

allo.Data <-
  within (allo.Data,
          {
            model.slopes <- summary (allo.mcmc.randomreg, random = TRUE) $ solutions
            model.slopes <- model.slopes [!grepl ('Intercept', rownames(model.slopes)), ]
            model.slopes <- model.slopes [-1, ]
            rownames (model.slopes) <-
              laply (strsplit(as.character (rownames (model.slopes)),
                              split = '\\.'), function (L) L [length (L)])

            model.slopes <-
              data.frame(model.slopes,
                         'node' = c (Tree $ NodeNames [-1, 2], names (OneDef)),
                         'pos' = c('Ancestor', 'Terminal') [
                           grepl ('_', rownames (model.slopes)) + 1])

            model.intercepts <- summary (allo.mcmc.randomreg, random = TRUE) $ solutions
            model.intercepts <- model.intercepts [!grepl ('logCS',
                                                          rownames(model.intercepts)), ]
            model.intercepts <- model.intercepts [-1, ]
            rownames (model.intercepts) <-
              laply (strsplit(as.character (rownames (model.intercepts)),
                              split = '\\.'), function (L) L [length (L)])

            model.intercepts <-
              data.frame(model.intercepts,
                         'node' = c (Tree $ NodeNames [-1, 2], names (OneDef)),
                         'pos' = c('Ancestor', 'Terminal') [
                           grepl ('_', rownames (model.intercepts)) + 1])

            
          })

allo.Plots $ bar.Wslope <-
  ggplot (allo.Data $ model.slopes) +
  geom_point(aes (x = node, y = post.mean, color = pMCMC < 0.05)) +
  geom_errorbar(aes (x = node, ymin = l.95..CI, ymax = u.95..CI), width = .1) +
  facet_wrap( ~ pos, ncol = 2, scales = 'free_x') +
### facet_grid(.id ~ pos, scales = 'free') +
  theme_bw() +
  theme(axis.text.x =
        element_text(angle = 90, size = 7, hjust = 1)) +
  geom_hline(xint = 0, linetype = 'dashed')

allo.Plots $ bar.Winter <-
  ggplot (allo.Data $ model.intercepts) +
  geom_point(aes (x = node, y = post.mean, color = pMCMC < 0.05)) +
  geom_errorbar(aes (x = node, ymin = l.95..CI, ymax = u.95..CI), width = .1) +
  facet_wrap( ~ pos, ncol = 2, scales = 'free_x') +
### facet_grid(.id ~ pos, scales = 'free') +
  theme_bw() +
  theme(axis.text.x =
        element_text(angle = 90, size = 7, hjust = 1)) +
  geom_hline(xint = 0, linetype = 'dashed')


allo.Data <-
  within (allo.Data,
          {
            phylo.slope.df <- model.slopes

            phylo.slope.df <-
              phylo.slope.df [match (c(Tree [[1]] $ tip.label, Tree $ NodeNames [-1, 2]),
                                    phylo.slope.df $ node),
                             c('node', 'post.mean', 'pMCMC', 'pos')]

            phylo.slope.df $ node <-
              factor (as.character (phylo.slope.df $ node),
                      levels = c(as.character (phylo.slope.df $ node) [c(108:216, 1:107)],
                        'Anthropoidea'))
            
            phylo.slope.df <- rbind (phylo.slope.df, c('Anthropoidea', 0, 1, 'Ancestor'))
            phylo.slope.df <- phylo.slope.df[c(1:109, 217, 110:216), ]

            phylo.slope.df $ post.mean <- as.numeric (phylo.slope.df $ post.mean)

            phylo.inter.df <- model.intercepts

            phylo.inter.df <-
              phylo.inter.df [match (c(Tree [[1]] $ tip.label, Tree $ NodeNames [-1, 2]),
                                    phylo.inter.df $ node),
                             c('node', 'post.mean', 'pMCMC', 'pos')]

            phylo.inter.df $ node <-
              factor (as.character (phylo.inter.df $ node),
                      levels = c(as.character (phylo.inter.df $ node) [c(108:216, 1:107)],
                        'Anthropoidea'))
            
            phylo.inter.df <- rbind (phylo.inter.df, c('Anthropoidea', 0, 1, 'Ancestor'))
            phylo.inter.df <- phylo.inter.df[c(1:109, 217, 110:216), ]

            phylo.inter.df $ post.mean <- as.numeric (phylo.inter.df $ post.mean)

          })

#### INTEGRATION

allo.Data <-
  within (allo.Data,
          {
            mod.index.local <-
              llply (W.all.ED, ModIndex, mod.hyp = Aux $ ed.hyp [[1]] [-20, 1:6],
                     it = 1000, .parallel = TRUE)
          })

allo.Data <-
  within (allo.Data,
          {
            mod.index.global <-
              llply (W.all.ED, ModIndex, mod.hyp = Aux $ ed.hyp [[1]] [-20, 7:8],
                     it = 1000, .parallel = TRUE)
          })

allo.Data <-
  within (allo.Data,
          {

            mod.index.df <- laply (mod.index.local,
                                   function (M) M [, 7])
            rownames (mod.index.df) <- names (W.all.Def)
            colnames (mod.index.df) <-
              c(paste (colnames (Aux $ def.hyp) [1:6], 'MI', sep = '.'), 'Full.MI')

            mod.index.df.gl <- laply (mod.index.global,
                                      function (M) M [, 7])
            rownames (mod.index.df.gl) <- names (W.all.Def)
            colnames (mod.index.df.gl) <-
              c(paste (colnames (Aux $ def.hyp) [7:8], 'MI', sep = '.'), 'Neuroface.MI')

            mod.index.df <- data.frame(mod.index.df)
            mod.index.df.gl <- data.frame(mod.index.df.gl)

            integra.df <-
              data.frame (
                laply (W.all.Def, ICV.Shape), ## ICV.Def, ICV.Def.Res
                laply (W.all.ED, ICV), ## ICV.ED, ICV.ED.Res
                c (laply (OneDef, function (L) sqrt (L $ ml.vcv [1, 1]) / L $ mean [1]),
                   laply (BW.Def, function (L)
                          sqrt (L $ W [1, 1]) /
                          weighted.mean (x = L $ means [, 1], w = L $ sample.sizes))), 
                mod.index.df, mod.index.df.gl
                )

            colnames (integra.df) [1:5] <- c('ICV.Def', 'ICV.Def.Res',
                                             'ICV.ED', 'ICV.ED.Res',
                                             'CSCV')
            
            rm (mod.index.df, mod.index.df.gl)

          })

allo.Data <-
  within (allo.Data,
          {

            wAC.all.Def <- laply (W.all.Def, function (C) c(C [-1, 1] / C [1, 1]))
            wAC.all.Def <- aaply (wAC.all.Def, 1, Normalize)
            
            rownames (wAC.all.Def) <- names (W.all.Def)
            colnames (wAC.all.Def) <- rownames (Aux $ def.hyp) [-1]
            
            integra.df $ wAC.CAC <- wAC.all.Def %*% root.CAC
            def.PC1 <- laply (W.all.Def, function (C) eigen (C[-1, -1]) $ vectors [, 1])

            integra.df $ wAC.PC1 <- rowSums (def.PC1 * wAC.all.Def)

            integra.df $ pos <-
              c ('Ancestor', 'Terminal') [grepl ('_', rownames (integra.df)) + 1]

            wAC.hyp <- wAC.all.Def %*% t (aaply (Aux $ def.hyp [-1, ], 2, Normalize))
            wAC.hyp <- data.frame (wAC.hyp)

            colnames (wAC.hyp) <- paste (colnames (Aux $ def.hyp), 'wAC', sep = '.')

            integra.df <- data.frame(integra.df, wAC.hyp,
                                     'W.slope' = phylo.slope.df $ post.mean,
                                     'W.inter' = phylo.inter.df $ post.mean)
            
          })

allo.Data $ integra.df $ node <- rownames (allo.Data $ integra.df)

allo.Data $ CAC.vs.hyp <- aaply (Aux $ def.hyp [-1, ], 2, Normalize) %*%
  allo.Data $ root.CAC

### MI Models

set.seed(23)

allo.Data $ A.inv <- inverseA(Tree [[1]]) $ Ainv

allo.Data <-
  within (allo.Data,
          {
            face.model <-
              MCMCglmm(Face.MI ~ W.slope,
                       random = ~ node,
                       rcov = ~ units, pr = TRUE, 
                       data = subset (integra.df, pos == 'Terminal'),
                       nitt = 70000, burnin = 40000, thin = 30,
                       ginverse = list (node = A.inv))
            
            neuro.model <-
              MCMCglmm(Neuro.MI ~ W.slope,
                       random = ~ node,
                       rcov = ~ units, pr = TRUE,
                       data = subset (integra.df, pos == 'Terminal'),
                       nitt = 70000, burnin = 40000, thin = 30,
                       ginverse = list (node = A.inv))

            oral.model <-
              MCMCglmm(Oral.MI ~ W.slope,
                       random = ~ node,
                       rcov = ~ units, pr = TRUE, 
                       data = subset (integra.df, pos == 'Terminal'),
                       nitt = 70000, burnin = 40000, thin = 30,
                       ginverse = list (node = A.inv))

            nasal.model <-
              MCMCglmm(Nasal.MI ~ W.slope,
                       random = ~ node,
                       rcov = ~ units, pr = TRUE,
                       data = subset (integra.df, pos == 'Terminal'),
                       nitt = 70000, burnin = 40000, thin = 30,
                       ginverse = list (node = A.inv))
            
            vault.model <-
              MCMCglmm(Vault.MI ~ W.slope,
                       random = ~ node,
                       rcov = ~ units, pr = TRUE,
                       data = subset (integra.df, pos == 'Terminal'),
                       nitt = 70000, burnin = 40000, thin = 30,
                       ginverse = list (node = A.inv))

            zygomatic.model <-
              MCMCglmm(Zygomatic.MI ~ W.slope,
                       random = ~ node,
                       rcov = ~ units, pr = TRUE,
                       data = subset (integra.df, pos == 'Terminal'),
                       nitt = 70000, burnin = 40000, thin = 30,
                       ginverse = list (node = A.inv))
          })

allo.Data <-
  within (allo.Data,
          {
            face.model.is <-
              MCMCglmm(Face.MI ~ W.slope + W.inter,
                       random = ~ node,
                       rcov = ~ units, pr = TRUE, 
                       data = subset (integra.df, pos == 'Terminal'),
                       nitt = 70000, burnin = 40000, thin = 30,
                       ginverse = list (node = A.inv))
            
            neuro.model.is <-
              MCMCglmm(Neuro.MI ~ W.slope + W.inter,
                       random = ~ node,
                       rcov = ~ units, pr = TRUE,
                       data = subset (integra.df, pos == 'Terminal'),
                       nitt = 70000, burnin = 40000, thin = 30,
                       ginverse = list (node = A.inv))

            oral.model.is <-
              MCMCglmm(Oral.MI ~ W.slope + W.inter,
                       random = ~ node,
                       rcov = ~ units, pr = TRUE, 
                       data = subset (integra.df, pos == 'Terminal'),
                       nitt = 70000, burnin = 40000, thin = 30,
                       ginverse = list (node = A.inv))

            nasal.model.is <-
              MCMCglmm(Nasal.MI ~ W.slope + W.inter,
                       random = ~ node,
                       rcov = ~ units, pr = TRUE,
                       data = subset (integra.df, pos == 'Terminal'),
                       nitt = 70000, burnin = 40000, thin = 30,
                       ginverse = list (node = A.inv))
            
            vault.model.is <-
              MCMCglmm(Vault.MI ~ W.slope + W.inter,
                       random = ~ node,
                       rcov = ~ units, pr = TRUE,
                       data = subset (integra.df, pos == 'Terminal'),
                       nitt = 70000, burnin = 40000, thin = 30,
                       ginverse = list (node = A.inv))

            zygomatic.model.is <-
              MCMCglmm(Zygomatic.MI ~ W.slope + W.inter,
                       random = ~ node,
                       rcov = ~ units, pr = TRUE,
                       data = subset (integra.df, pos == 'Terminal'),
                       nitt = 70000, burnin = 40000, thin = 30,
                       ginverse = list (node = A.inv))
          })

allo.Data <-
  within (allo.Data,
          {
            face.model.i <-
              MCMCglmm(Face.MI ~ W.inter,
                       random = ~ node,
                       rcov = ~ units, pr = TRUE, 
                       data = subset (integra.df, pos == 'Terminal'),
                       nitt = 70000, burnin = 40000, thin = 30,
                       ginverse = list (node = A.inv))
            
            neuro.model.i <-
              MCMCglmm(Neuro.MI ~ W.inter,
                       random = ~ node,
                       rcov = ~ units, pr = TRUE,
                       data = subset (integra.df, pos == 'Terminal'),
                       nitt = 70000, burnin = 40000, thin = 30,
                       ginverse = list (node = A.inv))

            oral.model.i <-
              MCMCglmm(Oral.MI ~ W.inter,
                       random = ~ node,
                       rcov = ~ units, pr = TRUE, 
                       data = subset (integra.df, pos == 'Terminal'),
                       nitt = 70000, burnin = 40000, thin = 30,
                       ginverse = list (node = A.inv))

            nasal.model.i <-
              MCMCglmm(Nasal.MI ~ W.inter,
                       random = ~ node,
                       rcov = ~ units, pr = TRUE,
                       data = subset (integra.df, pos == 'Terminal'),
                       nitt = 70000, burnin = 40000, thin = 30,
                       ginverse = list (node = A.inv))
            
            vault.model.i <-
              MCMCglmm(Vault.MI ~ W.inter,
                       random = ~ node,
                       rcov = ~ units, pr = TRUE,
                       data = subset (integra.df, pos == 'Terminal'),
                       nitt = 70000, burnin = 40000, thin = 30,
                       ginverse = list (node = A.inv))

            zygomatic.model.i <-
              MCMCglmm(Zygomatic.MI ~ W.inter,
                       random = ~ node,
                       rcov = ~ units, pr = TRUE,
                       data = subset (integra.df, pos == 'Terminal'),
                       nitt = 70000, burnin = 40000, thin = 30,
                       ginverse = list (node = A.inv))
          })



set.seed(42)

allo.Data <-
  within (allo.Data,
          {
            face.model.nophylo <-
              MCMCglmm(Face.MI ~ W.slope,
                       data = subset (integra.df, pos == 'Terminal'),
                       nitt = 70000, burnin = 40000, thin = 30)
                       
            
            neuro.model.nophylo <-
              MCMCglmm(Neuro.MI ~ W.slope,
                       data = subset (integra.df, pos == 'Terminal'),
                       nitt = 70000, burnin = 40000, thin = 30)

            oral.model.nophylo <-
              MCMCglmm(Oral.MI ~ W.slope,
                       data = subset (integra.df, pos == 'Terminal'),
                       nitt = 70000, burnin = 40000, thin = 30)

            nasal.model.nophylo <-
              MCMCglmm(Nasal.MI ~ W.slope,
                       data = subset (integra.df, pos == 'Terminal'),
                       nitt = 70000, burnin = 40000, thin = 30)
            
            vault.model.nophylo <-
              MCMCglmm(Vault.MI ~ W.slope,
                       data = subset (integra.df, pos == 'Terminal'),
                       nitt = 70000, burnin = 40000, thin = 30)
            
            zygomatic.model.nophylo <-
              MCMCglmm(Zygomatic.MI ~ W.slope,
                       data = subset (integra.df, pos == 'Terminal'),
                       nitt = 70000, burnin = 40000, thin = 30)
          })

### CI for regression

allo.Data <-
  within (allo.Data,
          {
            Ws.grid <- seq (min (subset (integra.df, pos == 'Terminal') $ W.slope),
                              max (subset (integra.df, pos == 'Terminal') $ W.slope),
                              length = 1000)
            
            face.pred <-
              aaply (face.model $ Sol, 1, function (ab) +
                     ab [1] + ab [2] * Ws.grid)
            
            face.pred <- adply(face.pred, 2, function (post)
                               c(mean (post),
                                 mean (post) - 1.96 * sd (post),
                                 mean (post) + 1.96 * sd (post)))

            face.pred $ Ws <- Ws.grid
            colnames (face.pred) <- c ('node', 'fit', 'lwr', 'upr', 'Ws')

            neuro.pred <-
              aaply (neuro.model $ Sol, 1, function (ab) +
                     ab [1] + ab [2] * Ws.grid)
            
            neuro.pred <- adply(neuro.pred, 2, function (post)
                                c(mean (post),
                                  mean (post) - 1.96 * sd (post),
                                  mean (post) + 1.96 * sd (post)))
            
            neuro.pred $ Ws <- Ws.grid
            colnames (neuro.pred) <- c ('node', 'fit', 'lwr', 'upr', 'Ws')

            oral.pred <-
              aaply (oral.model $ Sol, 1, function (ab) +
                     ab [1] + ab [2] * Ws.grid)
            
            oral.pred <- adply(oral.pred, 2, function (post)
                               c(mean (post),
                                 mean (post) - 1.96 * sd (post),
                                 mean (post) + 1.96 * sd (post)))
            
            oral.pred $ Ws <- Ws.grid
            colnames (oral.pred) <- c ('node', 'fit', 'lwr', 'upr', 'Ws')

            nasal.pred <-
              aaply (nasal.model $ Sol, 1, function (ab) +
                     ab [1] + ab [2] * Ws.grid)
            
            nasal.pred <- adply(nasal.pred, 2, function (post)
                                c(mean (post),
                                  mean (post) - 1.96 * sd (post),
                                  mean (post) + 1.96 * sd (post)))
            
            nasal.pred $ Ws <- Ws.grid
            colnames (nasal.pred) <- c ('node', 'fit', 'lwr', 'upr', 'Ws')

            vault.pred <-
              aaply (vault.model $ Sol, 1, function (ab) +
                     ab [1] + ab [2] * Ws.grid)
            
            vault.pred <- adply(vault.pred, 2, function (post)
                                c(mean (post),
                                  mean (post) - 1.96 * sd (post),
                                  mean (post) + 1.96 * sd (post)))
            
            vault.pred $ Ws <- Ws.grid
            colnames (vault.pred) <- c ('node', 'fit', 'lwr', 'upr', 'Ws')

            zygomatic.pred <-
              aaply (zygomatic.model $ Sol, 1, function (ab) +
                     ab [1] + ab [2] * Ws.grid)
            
            zygomatic.pred <- adply(zygomatic.pred, 2, function (post)
                                    c(mean (post),
                                      mean (post) - 1.96 * sd (post),
                                      mean (post) + 1.96 * sd (post)))
            
            zygomatic.pred $ Ws <- Ws.grid
            colnames (zygomatic.pred) <- c ('node', 'fit', 'lwr', 'upr', 'Ws')
            
          })


allo.Plots $ FaceMI.vs.Wslope <-
  ggplot (subset (allo.Data $ integra.df, pos == 'Terminal')) + 
  geom_point (aes (x = W.slope, y = Face.MI), shape = '+') +
  geom_smooth (aes (x = Ws, y = fit, ymin = lwr, ymax = upr),
               data = allo.Data $ face.pred, stat = 'identity', col = 'black') +
  theme_bw() + theme(text = element_text (size = 9)) +
  xlab (expression(b[s])) + ylab ('MI') + labs(title = 'Face')

allo.Plots $ NeuroMI.vs.Wslope <-
  ggplot (subset (allo.Data $ integra.df, pos == 'Terminal')) + 
  geom_point (aes (x = W.slope, y = Neuro.MI), shape = '+') +
  geom_smooth (aes (x = Ws, y = fit, ymin = lwr, ymax = upr),
               data = allo.Data $ neuro.pred, stat = 'identity', col = 'black') +
  theme_bw() + theme(text = element_text (size = 9)) +
  xlab (expression(b[s])) + ylab ('MI') + labs(title = 'Neuro')

  
allo.Plots $ OralMI.vs.Wslope <-
  ggplot (subset (allo.Data $ integra.df, pos == 'Terminal')) + 
  geom_point (aes (x = W.slope, y = Oral.MI), shape = '+') +
  geom_smooth (aes (x = Ws, y = fit, ymin = lwr, ymax = upr),
               data = allo.Data $ oral.pred, stat = 'identity', col = 'black') +
  theme_bw() + theme(text = element_text (size = 9)) +
  xlab (expression(b[s])) + ylab ('MI') + labs(title = 'Oral')

allo.Plots $ NasalMI.vs.Wslope <-
  ggplot (subset (allo.Data $ integra.df, pos == 'Terminal')) + 
  geom_point (aes (x = W.slope, y = Nasal.MI), shape = '+') +
  geom_smooth (aes (x = Ws, y = fit, ymin = lwr, ymax = upr),
               data = allo.Data $ nasal.pred, stat = 'identity', col = 'black') +
  theme_bw() + theme(text = element_text (size = 9)) +
  xlab (expression(b[s])) + ylab ('MI') + labs(title = 'Nasal')

allo.Plots $ VaultMI.vs.Wslope <-
  ggplot (subset (allo.Data $ integra.df, pos == 'Terminal')) + 
  geom_point (aes (x = W.slope, y = Vault.MI), shape = '+') +
  geom_smooth (aes (x = Ws, y = fit, ymin = lwr, ymax = upr),
               data = allo.Data $ vault.pred, stat = 'identity', col = 'black') +
  theme_bw() + theme(text = element_text (size = 9)) +
  xlab (expression(b[s])) + ylab ('MI') + labs(title = 'Vault')

allo.Plots $ ZygomaticMI.vs.Wslope <-
  ggplot (subset (allo.Data $ integra.df, pos == 'Terminal')) + 
  geom_point (aes (x = W.slope, y = Zygomatic.MI), shape = '+') +
  geom_smooth (aes (x = Ws, y = fit, ymin = lwr, ymax = upr),
               data = allo.Data $ zygomatic.pred, stat = 'identity', col = 'black') +
  theme_bw() + theme(text = element_text (size = 9)) +
  xlab (expression(b[s])) + ylab ('MI') + labs(title = 'Zygomatic')

allo.Plots $ MI.vs.Wslope <-
  arrangeGrob (arrangeGrob (allo.Plots $ FaceMI.vs.Wslope,
                            allo.Plots $ NeuroMI.vs.Wslope, ncol = 1),
               arrangeGrob (allo.Plots $ OralMI.vs.Wslope,
                            allo.Plots $ NasalMI.vs.Wslope,
                            allo.Plots $ ZygomaticMI.vs.Wslope,
                            allo.Plots $ VaultMI.vs.Wslope, ncol = 1),
               ncol = 2, widths = c(1, 0.6))

### LAMBDA!!!!

allo.Data $ lambda <- 
  cbind (allo.Data $ face.model $ VCV [, 1] /
         aaply (allo.Data $ face.model $ VCV, 1, sum),
         allo.Data $ neuro.model $ VCV [, 1] /
         aaply (allo.Data $ neuro.model $ VCV, 1, sum),
         allo.Data $ oral.model $ VCV [, 1] /
         aaply (allo.Data $ oral.model $ VCV, 1, sum),
         allo.Data $ nasal.model $ VCV [, 1] /
         aaply (allo.Data $ nasal.model $ VCV, 1, sum),
         allo.Data $ zygomatic.model $ VCV [, 1] /
         aaply (allo.Data $ zygomatic.model $ VCV, 1, sum),
         allo.Data $ vault.model $ VCV [, 1] /
         aaply (allo.Data $ vault.model $ VCV, 1, sum))

colnames (allo.Data $ lambda) <- c('Face', 'Neuro', 'Oral',
                                   'Nasal', 'Zygomatic', 'Vault')

allo.Data $ lambda <- melt (allo.Data $ lambda)

colnames (allo.Data $ lambda) <- c('it', 'hyp', 'lambda')

allo.Plots $ lambda.MI <- 
  ggplot (allo.Data $ lambda) +
  geom_boxplot (aes (y = lambda, x = hyp), outlier.shape = '+') +
  theme_bw() +
  xlab('Trait Set') +
  ylab (expression (lambda)) +
  theme(text = element_text (size = 9))
  
### source ('phylo.plot.R')


### 
attach ('../../Databases/Sym.RData')

### regress shape

allo.Data $ sym <- ldply (Sym, function (L) two.d.array(L $ sym))
allo.Data $ sym.3d <- arrayspecs (allo.Data $ sym [, -1], k = 3, p = 36)
allo.Data $ sym.gpa <- procGPA(allo.Data $ sym.3d, pcaoutput = FALSE)

allo.Data $ tan <- t (allo.Data $ sym.gpa $ tan)

allo.Data $ CAC.to.shape <-
  coef (lm (allo.Data $ tan ~ allo.Data $ onedef.df $ CACw)) [2, ]

Norm(allo.Data $ CAC.to.shape)

dim (allo.Data $ CAC.to.shape) <- c(36, 3)

allo.Data $ CAC.shapes <- array (0, c(36, 3, 5))
allo.Data $ cac.value <- c(-2, -1, 0, 1, 2)
for (i in 1:5)
  allo.Data $ CAC.shapes [, , i] <-
  allo.Data $ sym.gpa $ mshape + allo.Data $ cac.value [i] * allo.Data $ CAC.to.shape
  
dimnames (allo.Data $ CAC.shapes) [1:2] <- dimnames (Sym [[1]] $ sym) [1:2]

## source ('shape.plot.R')

## plot3d

for (i in 1:5)
  {
    coleurs <- colorRampPalette(c ('blue', 'red'))(5)
    points3d(allo.Data $ CAC.shape [, , i], rgl.open = ifelse (i == 1, T, F),
             col = coleurs [i])
    for (j in 1:(dim (Aux $ wireframe) [1]))
      lines3d (allo.Data $ CAC.shape [Aux $ wireframe [j, ], , i], col = coleurs [i])
  }

writeWebGL(dir = './', filename = 'CACshape.html', width = 768, height = 768)

browseURL('CACshape.html', 'firefox')

## tables

save (allo.Data, allo.Plots, file = 'allo.out.RData')

### CS vs ICV

allo.Data <-
  within (allo.Data,
          {
            size.model <-
              MCMCglmm(ICV.ED ~ CSCV,
                       random = ~ node,
                       rcov = ~ units, pr = TRUE,
                       data = subset (integra.df, pos == 'Terminal'),
                       nitt = 70000, burnin = 40000, thin = 30,
                       ginverse = list (node = A.inv))

            size.model.nophylo <-
              MCMCglmm(ICV.ED ~ CSCV,
                       data = subset (integra.df, pos == 'Terminal'),
                       nitt = 70000, burnin = 40000, thin = 30)
            
          })

summary (allo.Data $ size.model)
summary (allo.Data $ size.model.nophylo)


allo.Data <-
  within (allo.Data,
          {
            size.grid <- seq (min (subset (integra.df, pos == 'Terminal') $ CSCV),
                              max (subset (integra.df, pos == 'Terminal') $ CSCV),
                              length = 1000)

            size.pred <-
              aaply (size.model $ Sol, 1, function (ab) +
                     ab [1] + ab [2] * size.grid)

            size.pred <- adply(size.pred, 2, function (post)
                               c(mean (post),
                                 mean (post) - 1.96 * sd (post),
                                 mean (post) + 1.96 * sd (post)))

            size.pred $ sizevar <- size.grid
            colnames (size.pred) <- c ('node', 'fit', 'lwr', 'upr', 'sizevar')
          })

allo.Plots $ size.var <- 
  ggplot (subset (allo.Data $ integra.df, pos == 'Terminal')) + 
  geom_point (aes (x = CSCV, y = ICV.ED), shape = '+') +
  geom_smooth (aes (x = sizevar, y = fit, ymin = lwr, ymax = upr),
               data = allo.Data $ size.pred, stat = 'identity', col = 'black') +
  theme_bw() + theme(text = element_text (size = 9)) +
  xlab ('CSCV') + ylab('ICV')

allo.Data $ integra.dic <- 
  rbind (
    c (allo.Data $ face.model.nophylo $ DIC,
       allo.Data $ face.model $ DIC,
       allo.Data $ face.model.i $ DIC,
       allo.Data $ face.model.is $ DIC),
    c (allo.Data $ neuro.model.nophylo $ DIC,
       allo.Data $ neuro.model $ DIC,
       allo.Data $ neuro.model.i $ DIC,
       allo.Data $ neuro.model.is $ DIC),
    c (allo.Data $ oral.model.nophylo $ DIC,
       allo.Data $ oral.model $ DIC,
       allo.Data $ oral.model.i $ DIC,
       allo.Data $ oral.model.is $ DIC),
    c (allo.Data $ nasal.model.nophylo $ DIC,
       allo.Data $ nasal.model $ DIC,
       allo.Data $ nasal.model.i $ DIC,
       allo.Data $ nasal.model.is $ DIC),
    c (allo.Data $ zygomatic.model.nophylo $ DIC,
       allo.Data $ zygomatic.model $ DIC,
       allo.Data $ zygomatic.model.i $ DIC,
       allo.Data $ zygomatic.model.is $ DIC),
    c (allo.Data $ vault.model.nophylo $ DIC,
       allo.Data $ vault.model $ DIC,
       allo.Data $ vault.model.i $ DIC,
       allo.Data $ vault.model.is $ DIC)
    )

dimnames (allo.Data $ integra.dic) <- list (c('Face', 'Neuro', 'Oral',
                                              'Nasal', 'Zygomatic', 'Vault'),
                                            c('No Phylo', '$b_s$',
                                              '$a_s$', '$a_s + b_s$'))


allo.Data $ integra.dic <- allo.Data $ integra.dic [, -1]

## render ('allo.Rmd')
## system ('okular allo.pdf &')


## render ('supplemental.Rmd')
## system ('okular supplemental.pdf &')


save (allo.Data, allo.Plots, file = 'allo.Results.RData')

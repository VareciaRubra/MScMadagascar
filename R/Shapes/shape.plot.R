require (geomorph)
require (shapes)
require (RColorBrewer)
require (ape)
require (evolqg)
require (expm)
require (plyr)
require (plotrix)
require (doMC)
require (reshape2)
require (ggplot2)
require (phytools)
require (geiger)
require (mvtnorm)
require (MCMCglmm)
# require (grid)
# require (gridExtra)
# require (gridBase)
# require (surface)
require (rmarkdown)
require (knitr)
require (pander)
require (shape)
require (adephylo)
require (dplyr)
require (magrittr)
require (tidyr)
require (slidify)
require (slidifyLibraries)
require (StatMatch)
require (scales)
require (cowplot)
require (ggtree)
require (phylobase)

registerDoMC (cores = 3)

## require (devtools)
## install_github('muschellij2/slidify')
## install_github('ramnathv/slidifyLibraries')
## install_github('uyedaj/bayou')

attach ('~/MEGA/Doc/Databases/Reference.RData')
attach ('~/MEGA/Doc/Databases/ED.RData')
attach ('~/MEGA/Doc/Databases/Sym.RData')
attach ('~/MEGA/Doc/Databases/OneDef.RData')
attach ('~/MEGA/Doc/Databases/Tree.RData')
attach ('~/MEGA/Doc/Databases/Aux.RData')
attach ('~/MEGA/Doc/Databases/LifeHistory.RData')

options(contrasts = c('contr.sum', 'contr.poly'))

.source.files <- dir('~/MEGA/Doc/Tese/Func', pattern = '.R', full.names = TRUE)
.source.files <- .source.files [!grepl ('~', .source.files)]
for (i in 1:length (.source.files))
  source (.source.files [i])

attach('WPCs_Strepsirrhini.RData')
attach('WPCs_Strepsirrhini_ExtantsOnly.RData')

hapalemur <- read.csv('Hapalemur_7_AMNH_shape.csv', FALSE)

source ('~/MEGA/Doc/Func/glue.R')
source ('~/MEGA/Doc/Func/osymm.R')

glue.skull = function (Ai, Zi, sol = "svd")
      {
        mA = dim (Ai) [1]
        mZ = dim (Zi) [1]
        lmA = dimnames (Ai) [[1]]
        lmZ = dimnames (Zi) [[1]]
        cA = which (lmA %in% lmZ)
        cZ = which (lmZ %in% lmA)
        ### ordenar landmarks
        cA = cA[order (lmA[cA])]
        cZ = cZ[order (lmZ[cZ])]
        ### centroides dos pontos em comum
        ccA = t (array (colMeans (Ai[cA,]), c(3, mA)))
        ccZ = t (array (colMeans (Zi[cZ,]), c(3, mZ)))
        ### centralizando no centroide dos pontos em comum
        Ai = Ai - ccA
        Zi = Zi - ccZ
        ### angulos entre pontos em comum
        M = t (Ai[cA,]) %*% Zi[cZ,]
        UDV = svd (M)
        ### R é ortonormal, mas faz reflexões quando M faz
        D = diag (ifelse (UDV$d > 0, 1, -1))
        R = UDV$v %*% D %*% t(UDV$u)
        ### segunda rotação (sem expansão)
        Zi = Zi %*% R
        out = list (rbind (Ai, Zi [!(lmZ %in% lmA),]), det (R))
        return (out)
      }

dim (hapalemur)

hapa.lm <- hapalemur [, 1]
hapa.skull <- array (as.matrix (hapalemur [, 2:7]), c(27, 2, 3, 2))
hapa.skull <- aperm(hapa.skull, c(1, 3, 2, 4))

hapa.left <- hapa.skull [, , 1, ]
hapa.right <- hapa.skull [, , 2, ]

hapa.lm <- as.character (hapa.lm)
hapa.lm [3] <- hapa.lm [29] <- 'NA'
hapa.lm [22] <- 'TSe'

hapa.lm [28:54] <- gsub('e', 'd', hapa.lm[1:27])

dimnames (hapa.left) <- list (hapa.lm[1:27], c ('X', 'Y', 'Z'))
dimnames (hapa.right) <- list (hapa.lm[28:54], c ('X', 'Y', 'Z'))


hapa.full <- glue.skulls(hapa.left, hapa.right) [[1]]

hapa.av <- procGPA(hapa.full) $ mshape

dimnames (hapa.av) <- dimnames (hapa.full) [1:2]
rownames (hapa.full)

hapa.sym <- OSymm (hapa.av, 1:9, 28:45, 10:27) $ symmconf

rownames (hapa.sym) <- gsub ('d', '-D', gsub ('e', '-E', rownames (hapa.sym)))

hapa.sym <- hapa.sym [which (rownames (hapa.sym) %in% rownames (Reference)), ]

hapa.sym <- hapa.sym [match(rownames (Reference), rownames (hapa.sym)), ]

shape.plot <- list()

for (i in 1:4)
  shape.plot [[i]] <-
  ggshape(hapa.sym,
          Aux $ tessel.39 [1:39, ],
          WPCs[, i],
          rotation = c(1, -1, 1), culo = 0.01, thickness = 2) +
  ggtitle(i)

shape.plot [[5]] <- plot_grid(plotlist = shape.plot, ncol = 2)

shape.plot [[5]]

save (shape.plot, file = 'shape.plot.RData')

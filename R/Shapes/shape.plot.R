
load("~/MScMadagascar/R/Shapes/Aux.RData")
hapalemur <- read.csv('Hapalemur_7_AMNH_shape.csv', FALSE)

glue.skulls <- function (A, Z, soln = 'svd')
{
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
  if (dim (A) [3] != dim (Z) [3])
  {
    cat ('Número de vistas A não bate com número de vistas Z. Verifique!','\n')
    return (-1)
  }
  else
  {
    out = list ()
    dets = c()
    for (i in 1:dim (A)[3])
    {
      tmp = glue.skull (A[,,i], Z[,,i], sol = soln)
      if (i == 1)
        out = array (0, c(dim (tmp[[1]]), dim (A)[3]))
      out[,,i] = tmp[[1]]
      dets[i] = tmp[[2]] 
    }
    dimnames (out) = list (rownames (tmp[[1]]),
                           colnames (tmp[[1]]),
                           dimnames (A) [[3]])
    right.first = function (element)
    {
      return (ifelse (length (element) == 1, 0,
                      ifelse (element [1] == 'NLT', 3,
                              ifelse (element [2] == 'D', 1, 2))))
    }
    ord = strsplit (dimnames (out) [[1]], split = '-')
    ord = sapply (ord, right.first)
    out = out [order (ord),,]
    return (list (out, dets))
  }
}

Reference <- c("IS","NSL","NA", "BR", "PNS","BA", "OPI","LD", "PM-D","PT-D","FM-D","ZS-D","ZI-D","MT-D","APET-D","EAM-D","PEAM-D","ZYGO-D","TSP-D","TS-D","JP-D","AS-D","PM-E","PT-E","FM-E","ZS-E","ZI-E","MT-E","APET-E","EAM-E","PEAM-E","ZYGO-E", "TSP-E","TS-E","JP-E","AS-E")
SHAPE = hapalemur
W.MATRIX = W.matrix

PCLoadShapePlotter <- function( SHAPE, W.MATRIX )
{

  dim (SHAPE)
  
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


WPCs <- eigen(W.MATRIX)$vectors
WPCs[, 1] <- -WPCs[, 1]

traits <- rownames(SHAPE.sym) %>% strsplit(split = "-") %>% laply(function (x) x[1])

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
shape.plot <- list()
for (i in 1:39)
{
myBreaks <- c(round(min(WPCs[, i]), digits = 1), round(mean (WPCs[, i]), digits = 1 ), round(max(WPCs[, i]), digits = 1 )) 
interval <- mean (WPCs[, i])
shape.plot [[i]] <-
  ggshape(shape = SHAPE.sym,
          wireframe = Aux $ tessel.39 [1:39, ],
          colors = WPCs[, i],
          rotation = c(1, -1, 1), 
          culo = 0.02, 
          thickness = 2) +
  geom_point (aes (x = X, y = Y), alpha = 0.6) +
  geom_label(aes (x = X, y = Y, label = traits ),  alpha = 0.6) +
  ggtitle(paste("PC", i, sep = " ")) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position= c(0.3,0.1), legend.direction="horizontal") +
  scale_fill_gradientn("", limits = c(min(WPCs[, i])-0.09, max(WPCs[, i])+0.09),
                        breaks= myBreaks,
                        guide = guide_colorbar(nbin=100, draw.ulim = T, draw.llim = T),
                        colors = myPalette(11)) +
  scale_color_gradientn("", limits = c(min(WPCs[, i]) -0.09, max(WPCs[, i]) +0.09 ),
                        breaks= myBreaks,
                        guide = guide_colorbar(nbin=100, draw.ulim = T, draw.llim = T),
                        colors = myPalette(11))
}
shape.plot [[4]]
shape.plot [[40]] <- plot_grid(plotlist = shape.plot, ncol = 4) + ggtitle("Trais scores on W\nusing Hapalemur frame")
shape.plot [[41]] <- plot_grid(shape.plot[[1]], shape.plot[[2]],shape.plot[[3]], shape.plot[[4]], ncol = 2) + ggtitle("Trais scores on W\nusing Hapalemur frame")
shape.plot [[41]]
#rm(shape.plot)

return( list("W.Loadings.Shape" = shape.plot ))

}

temp <- PCLoadShapePlotter(SHAPE = hapalemur, W.MATRIX = W.matrix)
temp$W.Loadings.Shape[[41]]

#save (shape.plot, file = 'shape.plot.RData')

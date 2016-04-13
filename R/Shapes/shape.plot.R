
load("~/MScMadagascar/R/Shapes/Aux.RData")
tupaia <- read.csv('R/Shapes/shapes/Tupaia_4A_RMNH_shape.csv', FALSE)
hapalemur <- read.csv('R/Shapes/shapes/Hapalemur_7_AMNH_shape.csv', FALSE)
varecia <- read.csv('R/Shapes/shapes/Varecia_8B_RMNH_shape.csv', FALSE)
microcebus <- read.csv('R/Shapes/shapes/Microcebus_1A_MNHN_shape.csv', FALSE)
lepilemur <- read.csv('R/Shapes/shapes/Lepilemur_1A_RMNH_shape.csv', FALSE)
propithecus <- read.csv('R/Shapes/shapes/Propithecus_1A_RMNH_shape.csv', FALSE)
daubentonia <- read.csv('R/Shapes/shapes/Daubentonia_4A_RMNH_shape.csv', FALSE)



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
c(1, -1, 1)

PCLoadShapePlotter <- function( SHAPE, W.MATRIX, ROTACIONI = c(1, -1, 1))
  
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
          rotation = ROTACIONI , 
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

PC.Plots <- vector("list")

PC.Plots$Tarsius <- ancestral.mx$'43'
PC.Plots$Microcebus <- PCLoadShapePlotter(SHAPE = microcebus, W.MATRIX = mx.list.taxonomy$Microcebus, ROTACIONI = c(1,-1,1))
PC.Plots$Mirza <- ancestral.mx$Mirza_coquereli
PC.Plots$Cheirogaleus <- ancestral.mx$'52'
PC.Plots$Phaner <- ancestral.mx$Phaner_furcifer
PC.Plots$W.Cheirogaleidae <- ancestral.mx$'48'
PC.Plots$Lepilemur <- PCLoadShapePlotter(SHAPE = lepilemur, W.MATRIX = mx.list.taxonomy$Lepilemur, ROTACIONI = c(-1,-1,1))
PC.Plots$Avahi <- ancestral.mx$'58'
PC.Plots$Propithecus <- PCLoadShapePlotter(SHAPE = propithecus, W.MATRIX = mx.list.taxonomy$Propithecus, ROTACIONI = c(-1,-1,1))
PC.Plots$Indri <- ancestral.mx$Indri_indri
PC.Plots$W.Indridae <- ancestral.mx$'56'
PC.Plots$Eulemur <- ancestral.mx$'65'
PC.Plots$Hapalemur <- ancestral.mx$Hapalemur_griseus
PC.Plots$Prolemur <- ancestral.mx$Prolemur_simus
PC.Plots$Lemur <- ancestral.mx$Lemur_catta
PC.Plots$Varecia <- PCLoadShapePlotter(SHAPE = varecia, W.MATRIX = mx.list.taxonomy$Varecia, ROTACIONI = c(-1,-1,1))
PC.Plots$W.Lemuridae <- ancestral.mx$'63'
PC.Plots$Daubentonia <- PCLoadShapePlotter(SHAPE = daubentonia, W.MATRIX = mx.list.taxonomy$Daubentonia, ROTACIONI = c(1,-1,1))
PC.Plots$W.Madagascar <- ancestral.mx$'45'
PC.Plots$Perodicticus <- ancestral.mx$Perodicticus_potto
PC.Plots$Loris <- ancestral.mx$Loris_tardigradus
PC.Plots$W.Lorisidae <- ancestral.mx$'78'
PC.Plots$Nycticebus <- ancestral.mx$Nycticebus_coucang
PC.Plots$Euoticus <- ancestral.mx$Euoticus_elegantulus
PC.Plots$Otolemur <- ancestral.mx$Otolemur_crassicaudatus
PC.Plots$Galago <- ancestral.mx$Galago_senegalensis
PC.Plots$W.Galagidae <- ancestral.mx$'80'
PC.Plots$W.OutMadagascar <- ancestral.mx$'77'
PC.Plots$W.Strepsirrhini <-PCLoadShapePlotter(SHAPE = hapalemur, W.MATRIX = mx.list.taxonomy$W.Strepsirrhini, ROTACIONI = c(1,-1,1))
PC.Plots$W.Prosimian <-PCLoadShapePlotter(SHAPE = tupaia, W.MATRIX = mx.list.taxonomy$W.Prosimian, ROTACIONI = c(-1,-1,1))


PC.Plots <- 
temp$W.Loadings.Shape[[41]]

save (PC.Plots, file = 'mx_list_taxonomy.RData')

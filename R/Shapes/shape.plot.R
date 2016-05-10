load("~/MScMadagascar/R/Shapes/Aux.RData")
load("~/MScMadagascar/R/Shapes/shapes/mx_list_taxonomy.RData")
load("~/MScMadagascar/R/Shapes/file.RData")
load("~/MScMadagascar/R/Shapes/shapes/traits.RData")

source('~/MScMadagascar/R/LoadPacks.R', echo=TRUE)


tupaia <- read.csv('R/Shapes/shapes/Tupaia_4A_RMNH_shape.csv', FALSE)
allocebus <- read.csv('R/Shapes/shapes/Allocebus_1A_MNHN_shape.csv', FALSE)
microcebus <- read.csv('R/Shapes/shapes/Microcebus_1A_MNHN_shape.csv', FALSE)
cheirogaleus <- read.csv('R/Shapes/shapes/Cheirogaleus_2B_MNHN_shape.csv', FALSE)
phaner <- read.csv('R/Shapes/shapes/Phaner_6A_RMNH_shape.csv', FALSE)
mirza <- read.csv('R/Shapes/shapes/Mirza_1A_RMNH_shape.csv', FALSE)
varecia <- read.csv('R/Shapes/shapes/Varecia_8B_RMNH_shape.csv', FALSE)
eulemur <- read.csv('R/Shapes/shapes/Eulemur_7G_AMNH_shape.csv', FALSE)
prolemur <- read.csv("R/Shapes/shapes/Prolemur_1A_ZMB_shape.csv", FALSE)
hapalemur <- read.csv('R/Shapes/shapes/Hapalemur_7_AMNH_shape.csv', FALSE)
lemur <- read.csv('R/Shapes/shapes/Lemur_2A_AMNH_shape.csv', FALSE)
lepilemur <- read.csv('R/Shapes/shapes/Lepilemur_8A_AMNH_shape.csv', FALSE)
indri <- read.csv('R/Shapes/shapes/Indri_5A_RMNH_shape.csv', FALSE)
propithecus <- read.csv('R/Shapes/shapes/Propithecus_1A_RMNH_shape.csv', FALSE)
avahi <- read.csv('R/Shapes/shapes/Avahi_10A_MNHN_shape.csv', FALSE)
daubentonia <- read.csv('R/Shapes/shapes/Daubentonia_4A_RMNH_shape.csv', FALSE)
euoticus <- read.csv('R/Shapes/shapes/Euoticus_2A_AMNH_shape.csv', FALSE)
otolemur <- read.csv('R/Shapes/shapes/Otolemur_1C_USNM_shape.csv', FALSE)
galago <- read.csv('R/Shapes/shapes/Galago_1C_USNM_shape.csv', FALSE)
loris <- read.csv('R/Shapes/shapes/Loris_1A_FMNH_shape.csv', FALSE)
nycticebus <- read.csv('R/Shapes/shapes/Nycticebus_3B_FMNH_shape.csv', FALSE)
perodicticus <- read.csv('R/Shapes/shapes/Perodicticus_10A_FMNH_shape.csv', FALSE)
arctocebus <- read.csv('R/Shapes/shapes/Arctocebus_1A_FMNH_shape.csv', FALSE)
tarsius <- read.csv('R/Shapes/shapes/Tarsius_1B_FMNH_shape.csv', FALSE)

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

PCLoadShapePlotter <- function( SHAPE, W.MATRIX, ROTACIONI = c(1, -1, 1), TTL = "Awesome Lemur")
  
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
if (WPCs[[1]] <=0) WPCs <- -WPCs
if (WPCs[[1]] >=0) WPCs <- WPCs

#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
myPalette <- colorRampPalette(c("red", "grey", "blue") )
shape.plot <- list()
for (i in 1:6)
{
myBreaks <- c(round(min(WPCs[, i]) , digits = 1), round(mean (WPCs[, i]), digits = 1 ), round(max(WPCs[, i]), digits = 1 )) 
interval <- mean (WPCs[, i])
shape.plot [[i]] <-
  ggshape(shape = SHAPE.sym,
          wireframe = Aux $ tessel.39 [1:39, ],
          colors = rev(WPCs[, i]),
          rotation = ROTACIONI , 
          culo = 0.02, 
          thickness = 0.9) +
  #geom_point (aes (x = X, y = Y), alpha = 0.6, color ="darkgrey", size = 1.2) +
  #geom_label(aes (x = X, y = Y, label = traits ),  alpha = 0.6, size = 7, label.padding = unit(0.4, "mm")) +
  #ggtitle(paste("PC", i, sep = " ")) +
  theme(plot.title = element_text(face = "bold", size = 12),
        #legend.position = "none",
        legend.position= c(0.3,0.1), 
        #legend.direction="horizontal",
        legend.text = element_text(size = 8),
        plot.margin = unit(c(0,0,-0,0), "cm"), 
        legend.key.size = unit(0.3, "cm"), 
        panel.margin.x = unit(0.3, "cm"), 
        panel.margin.y = unit(0.3, "cm") ) +
  scale_fill_gradient2(paste("PC", i, sep = " "),  low = "red", mid = "grey", high ="blue",  midpoint = 0, 
                       #labels =c(-0.7, 0, 0.7), 
                       #breaks = c(-0.7, 0, 0.7),
                        breaks= myBreaks,
                        #colors = myPalette(11),
                        guide = guide_colorbar(nbin=102, draw.ulim = T, draw.llim = T) ) +
  scale_color_gradient2(paste("PC", i, sep = " "), low = "red", mid = "grey", high ="blue",  midpoint = 0,
                        breaks= myBreaks,
                        #colors = myPalette(11),
                        guide = guide_colorbar(nbin=102, draw.ulim = T, draw.llim = T) ) 
rm(myBreaks)}
#shape.plot [[40]] <- plot_grid(plotlist = shape.plot, ncol = 4, labels = TTL) 
shape.plot [[7]] <- plot_grid(shape.plot[[1]], shape.plot[[2]],shape.plot[[3]],shape.plot[[4]], ncol =  4, align = "hv", labels = TTL) 
shape.plot [[8]] <- shape.plot[[1]] + scale_fill_gradientn(colours = "darkgrey") + scale_color_gradientn(colours = "darkgrey") + ggtitle(TTL) + theme(plot.title  = element_text(face = "italic", size = 15), legend.position = "none" )
#rm(shape.plot)
return( shape.plot)
}

temp <- PCLoadShapePlotter(SHAPE = hapalemur, W.MATRIX = mx.list.taxonomy$W.Madagascar, ROTACIONI = c(1,-1,1), TTL = "Strepsirrhini W matrix")
temp[[7]]

PC.Plots$Galago <-PCLoadShapePlotter(SHAPE = galago, W.MATRIX = mx.list.taxonomy$Galago, ROTACIONI = c(-1,-1,1), TTL = "Galago")
PC.Plots$Galago [[8]]
PC.Plots <- vector("list")

PC.Plots$Galago <- PCLoadShapePlotter(SHAPE = galago, W.MATRIX = mx.list.taxonomy$Galago, ROTACIONI = c(-1,-1,1), TTL = "Galago")
PC.Plots$Otolemur <- PCLoadShapePlotter(SHAPE = otolemur, W.MATRIX = mx.list.taxonomy$Otolemur, ROTACIONI = c(1,-1,1), TTL = "Otolemur")
PC.Plots$Euoticus <- PCLoadShapePlotter(SHAPE = euoticus, W.MATRIX = mx.list.taxonomy$Euoticus, ROTACIONI = c(-1,-1,1), TTL = "Euoticus")
PC.Plots$Nycticebus <- PCLoadShapePlotter(SHAPE = nycticebus, W.MATRIX = mx.list.taxonomy$W.Lorisidae, ROTACIONI = c(-1,-1,1), TTL = "Nycticebus")
PC.Plots$Loris <- PCLoadShapePlotter(SHAPE = loris, W.MATRIX = mx.list.taxonomy$W.Lorisidae, ROTACIONI = c(-1,-1,1), TTL = "Loris")
PC.Plots$Perodicticus <- PCLoadShapePlotter(SHAPE = perodicticus, W.MATRIX = mx.list.taxonomy$Perodicticus, ROTACIONI = c(-1,-1,1), TTL = "Perodicticus")
PC.Plots$Arctocebus <- PCLoadShapePlotter(SHAPE = arctocebus, W.MATRIX = mx.list.taxonomy$W.Lorisidae, ROTACIONI = c(1,-1,1), TTL = "Arctocebus")

PC.Plots$Daubentonia <- PCLoadShapePlotter(SHAPE = daubentonia, W.MATRIX = mx.list.taxonomy$Daubentonia, ROTACIONI = c(1,-1,1), TTL = "Daubentonia")

PC.Plots$Varecia <- PCLoadShapePlotter(SHAPE = varecia, W.MATRIX = mx.list.taxonomy$Varecia, ROTACIONI = c(-1,-1,1), TTL = "Varecia")
PC.Plots$Lemur <- PCLoadShapePlotter(SHAPE = lemur, W.MATRIX = mx.list.taxonomy$Lemur, ROTACIONI = c(1,-1,1), TTL = "Lemur")
PC.Plots$Prolemur <- PCLoadShapePlotter(SHAPE = prolemur, W.MATRIX = mx.list.taxonomy$Prolemur, ROTACIONI = c(-1,-1,1), TTL = "Prolemur")
PC.Plots$Hapalemur <- PCLoadShapePlotter(SHAPE = hapalemur, W.MATRIX = mx.list.taxonomy$Hapalemur, ROTACIONI = c(1,-1,1), TTL = "Hapalemur")
PC.Plots$Eulemur <- PCLoadShapePlotter(SHAPE = eulemur, W.MATRIX = mx.list.taxonomy$Eulemur, ROTACIONI = c(1,-1,1), TTL = "Eulemur")


PC.Plots$Indri <- PCLoadShapePlotter(SHAPE = indri, W.MATRIX = mx.list.taxonomy$Indri, ROTACIONI = c(-1,-1,1), TTL = "Indri")
PC.Plots$Propithecus <- PCLoadShapePlotter(SHAPE = propithecus, W.MATRIX = mx.list.taxonomy$Propithecus, ROTACIONI = c(-1,-1,1), TTL = "Propithecus")
PC.Plots$Avahi <- PCLoadShapePlotter(SHAPE = avahi, W.MATRIX = mx.list.taxonomy$Avahi_laniger, ROTACIONI = c(-1,-1,1), TTL ="Avahi")

PC.Plots$Lepilemur <- PCLoadShapePlotter(SHAPE = lepilemur, W.MATRIX = mx.list.taxonomy$Lepilemur, ROTACIONI = c(1,-1,1), TTL = "Lepilemur")
PC.Plots$Phaner <- PCLoadShapePlotter(SHAPE = phaner, W.MATRIX = mx.list.taxonomy$Phaner, ROTACIONI = c(-1,-1,1), TTL ="Phaner")
PC.Plots$Cheirogaleus <- PCLoadShapePlotter(SHAPE = cheirogaleus, W.MATRIX = mx.list.taxonomy$Cheirogaleus, ROTACIONI = c(1,-1,1), TTL = "Cheirogaleus")
PC.Plots$Mirza <- PCLoadShapePlotter(SHAPE = mirza, W.MATRIX = mx.list.taxonomy$Mirza, ROTACIONI = c(-1,-1,1), TTL = "Mirza")
PC.Plots$Microcebus <- PCLoadShapePlotter(SHAPE = microcebus, W.MATRIX = mx.list.taxonomy$Microcebus, ROTACIONI = c(1,-1,1), TTL = "Microcebus")
PC.Plots$Allocebus<- PCLoadShapePlotter(SHAPE = allocebus, W.MATRIX = mx.list.taxonomy$Microcebus, ROTACIONI = c(1,-1,1), TTL = "Allocebus")
PC.Plots$Tarsius <- PCLoadShapePlotter(SHAPE = tarsius, W.MATRIX = mx.list.taxonomy$Tarsius, ROTACIONI = c(-1,-1,1), TTL = "Tarsius")

PC.Plots$W.Daubentonia <- PCLoadShapePlotter(SHAPE = daubentonia, W.MATRIX = mx.list.taxonomy$Daubentonia, ROTACIONI = c(1,-1,1), TTL = "Daubentonia (g)")
PC.Plots$W.Lemuridae <- PCLoadShapePlotter(SHAPE = lemur, W.MATRIX = mx.list.taxonomy$W.Lemuridae, ROTACIONI = c(1,-1,1), TTL ="Lemuridae (f)")
PC.Plots$W.Indridae <- PCLoadShapePlotter(SHAPE = propithecus, W.MATRIX = mx.list.taxonomy$W.Indridae, ROTACIONI = c(-1,-1,1), TTL = "Indridae (f)")
PC.Plots$W.Lepilemuridae <- PCLoadShapePlotter(SHAPE = lepilemur, W.MATRIX = mx.list.taxonomy$Lepilemur, ROTACIONI = c(1,-1,1), TTL = "Lepilemuridae (f)")
PC.Plots$W.Cheirogaleidae <- PCLoadShapePlotter(SHAPE = microcebus, W.MATRIX = mx.list.taxonomy$W.Cheirogaleidae, ROTACIONI = c(1,-1,1), TTL ="Cheirogaleidae (f)")
PC.Plots$W.Madagascar <- PCLoadShapePlotter(SHAPE = hapalemur, W.MATRIX = mx.list.taxonomy$W.Madagascar, ROTACIONI = c(1,-1,1), TTL = "Malagasy clade")

PC.Plots$W.Lorisidae <- PCLoadShapePlotter(SHAPE = loris, W.MATRIX = mx.list.taxonomy$W.Lorisidae, ROTACIONI = c(-1,-1,1), TTL = "Lorisidae (f)")
PC.Plots$W.Galagidae <- PCLoadShapePlotter(SHAPE = galago, W.MATRIX = mx.list.taxonomy$W.Galagidae, ROTACIONI = c(-1,-1,1), TTL = "Galagidae (f)")
PC.Plots$W.OutMadagascar <- PCLoadShapePlotter(SHAPE = galago, W.MATRIX = mx.list.taxonomy$W.OutMadagascar, ROTACIONI = c(-1,-1,1), TTL = "Lorisiformes")
PC.Plots$W.Strepsirrhini <- PCLoadShapePlotter(SHAPE = hapalemur, W.MATRIX = mx.list.taxonomy$W.Strepsirrhini, ROTACIONI = c(1,-1,1), TTL = "Strepsirrhini")
PC.Plots$Tarsius.g <- PCLoadShapePlotter(SHAPE = tarsius, W.MATRIX = mx.list.taxonomy$Tarsius, ROTACIONI = c(-1,-1,1), TTL = "Tarsius (g)")
PC.Plots$W.Prosimian <-PCLoadShapePlotter(SHAPE = tupaia, W.MATRIX = mx.list.taxonomy$W.Prosimian, ROTACIONI = c(-1,-1,1), TTL = "Tupaia (g)")

save.image(file = 'Data/PC_loadings_plots.RData')

shapes.general <- PC.Plots[1:7] %>% llply(function (x) x[[8]]) 
plot.shapes.out.madagascar <- plot_grid(plotlist = shapes.general, ncol = 1)
save_plot(filename = "R/Shapes/plot.shapes.out.madagascar.pdf", plot = plot.shapes.out.madagascar, 
          base_aspect_ratio = 0.9, base_height = 14, base_width = 4)

shapes.general <- PC.Plots[8:13] %>% llply(function (x) x[[8]]) 
plot.shapes.lemuridae<- plot_grid(plotlist = shapes.general, ncol = 1)
save_plot(filename = "R/Shapes/plot.shapes.lemuridae.pdf", plot = plot.shapes.lemuridae, 
          base_aspect_ratio = 0.9, base_height = 14, base_width = 4)

shapes.general <- PC.Plots[14:17] %>% llply(function (x) x[[8]]) 
plot.shapes.indridae <- plot_grid(plotlist = shapes.general, ncol = 1)
save_plot(filename = "R/Shapes/plot.shapes.indridae.pdf", plot = plot.shapes.indridae, 
          base_aspect_ratio = 0.9, base_height = 14, base_width = 4)

shapes.general <- PC.Plots[18:23] %>% llply(function (x) x[[8]]) 
plot.shapes.cheirohgaleidae <- plot_grid(plotlist = shapes.general, ncol = 1)
save_plot(filename = "R/Shapes/plot.shapes.cheirogaleidae.pdf", plot = plot.shapes.cheirohgaleidae, 
          base_aspect_ratio = 0.9, base_height = 14, base_width = 4)

shapes.general <- PC.Plots[24:28] %>% llply(function (x) x[[8]]) 
plot.shapes.madagascar <- plot_grid(plotlist = shapes.general, ncol = 1)
save_plot(filename = "R/Shapes/plot.shapes.madagascar.pdf", plot = plot.shapes.madagascar, 
          base_aspect_ratio = 0.9, base_height = 14, base_width = 4)

shapes.general <- PC.Plots[30:35] %>% llply(function (x) x[[8]]) 
plot.shapes.prosimian <- plot_grid(plotlist = shapes.general, ncol = 1)
save_plot(filename = "R/Shapes/plot.shapes.prosimian.pdf", plot = plot.shapes.prosimian, 
          base_aspect_ratio = 0.9, base_height = 14, base_width = 4)

shapes.general <- PC.Plots[1:7] %>% llply(function (x) x[[7]]) 
plot.pc.out.madagascar <- plot_grid(plotlist = shapes.general, ncol = 1)
save_plot(filename = "R/Shapes/plot.pc.out.madagascar.pdf", plot = plot.pc.out.madagascar, 
          base_aspect_ratio = 0.9, base_height = 14, base_width = 7)


shapes.PC.main.Madagascar <-plot_grid(PC.Plots$Daubentonia[[7]], PC.Plots$W.Lemuridae[[7]], PC.Plots$W.Indridae[[7]], PC.Plots$W.Lepilemuridae[[7]], PC.Plots$W.Cheirogaleidae[[7]], ncol = 1, align = "hv")
save_plot(filename = "R/Shapes/shapes.PC.main.Madagascar.pdf", plot = shapes.PC.main.Madagascar, 
          base_aspect_ratio = 1, base_height = 12)

shapes.PC.main.Out.Madagascar <-plot_grid(PC.Plots$W.Galagidae[[7]], PC.Plots$W.Lorisidae[[7]], PC.Plots$W.OutMadagascar[[7]], ncol = 1, align = "hv")
save_plot(filename = "R/Shapes/shapes.PC.main.OutMadagascar.pdf", plot = shapes.PC.main.Out.Madagascar, 
          base_aspect_ratio = 1, base_height = 12)



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

Shapes <- vector("list")
Shapes$tupaia <- tupaia 
Shapes$allocebus <- allocebus
Shapes$microcebus <- microcebus
Shapes$cheirogaleus <- cheirogaleus
Shapes$phaner <- phaner
Shapes$mirza <- mirza
Shapes$varecia <- varecia
Shapes$eulemur <- eulemur
Shapes$prolemur <- prolemur
Shapes$hapalemur <- hapalemur
Shapes$lemur <- lemur
Shapes$lepilemur <- lepilemur
Shapes$indri <- indri
Shapes$propithecus <- propithecus
Shapes$avahi <- avahi
Shapes$daubentonia <- daubentonia
Shapes$euoticus <- euoticus
Shapes$otolemur <- otolemur
Shapes$galago <- galago
Shapes$loris <- loris
Shapes$nycticebus <- nycticebus
Shapes$perodicticus <- perodicticus
Shapes$arctocebus <- arctocebus
Shapes$tarsius <- tarsius


Shapes.sym <- Shapes %>% llply(Shaper)
save(Shapes.sym, file = "Data/Shapes.RData")

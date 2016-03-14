ggshape <- 
function (shape, wireframe, colors, view = c(1, 2, 3),
          rotation = c(1, 1, 1), culo = 0.015, thickness = 3,
          palette = rev (brewer.pal (10, 'Spectral')))
{
  Q <- shape
  Lms <- rownames (Q)
  Right <- grep ('-D', Lms)
  Left <- grep ('-E', Lms)
  Midline <- !(Lms %in% c (Lms [Right], Lms [Left]))
  Q <- Q %*% Rotate2MidlineMatrix (Q, Lms [Midline]) 
  Q <- Q %*% diag (rotation)
  if (all (view == c(1, 2, 3)))
    Q <- Q %*%
    array (c(cos(pi/7), sin (pi/7), 0,
             -sin (pi/7), cos(pi/7), 0,
             0, 0, 1), c(3, 3))
  colnames (Q) <- c('X', 'Y', 'Z')
  ## pts <- which (rownames (Q) %in% rownames (Q.tetra))
  Q <- Q [, view]
  colnames (Q) <- c('X', 'Y', 'Z')
  Q.tetra <- Q [wireframe, ]
  dim (Q.tetra) <- c (dim (wireframe), ncol (Q))
  Q.names <- array (rownames (Q) [wireframe], dim = dim (wireframe))
  Q.names <- apply (Q.names, 1, paste, collapse = '.')
  dimnames (Q.tetra) <- list('ild' = Q.names,
                             'pos' = c(1, 2),
                             'dim' = c('X', 'Y', 'Z'))
  
  Q.singular <- which (!duplicated (dimnames (Q.tetra) [[1]]))
  
  Q.colors <- rep (colors, times = 2) [Q.singular]
  
  Q.tetra <- Q.tetra [Q.singular, , ]
  
  Q.tetra.df <- dcast (melt (Q.tetra), ild ~ dim + pos)
  
  Q.tetra.df $ culo <- rep (culo, nrow (Q.tetra))
  
  #Q.tetra.df $ Z_1 <-
  #  (Q.tetra.df $ Z_1 - min (Q.tetra.df $ Z_1)) /
  #    (max (Q.tetra.df $ Z_1) - min (Q.tetra.df $ Z_1))
  
  #Q.tetra.df $ Z_2 <-
  #  (Q.tetra.df $ Z_2 - min (Q.tetra.df $ Z_2)) /
  #    (max (Q.tetra.df $ Z_2) - min (Q.tetra.df $ Z_2))
  
  Q.line.df <-
    plyr::ddply (Q.tetra.df, .(ild), plyr::summarise,
                 'X' = c(
                   X_1 - Z_1 * culo,
                   X_1 + Z_1 * culo,
                   X_2 + Z_2 * culo,
                   X_2 - Z_2 * culo),
                 'Y' = c(
                   Y_1,
                   Y_1,
                   Y_2,
                   Y_2))
  
  Q.line.df $ color <- rep (Q.colors, each = 4)
  
  shape.plot <-
    ggplot (data.frame(Q)) +
    geom_point (aes (x = X, y = Y), alpha = 0) +
    coord_fixed() +
    theme_minimal() +
    theme(plot.margin = unit(c(0, 0, 0, 0), 'cm')) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    guides(size = FALSE) +
    xlab('') + ylab('')
  
  nominhos <- rownames(Q)
  nominhos <- nominhos %>% gsub("-D", "", .) %>% gsub("-E", "",.)
  
  breaks.values <- as.vector(round(summary(colors)[c(1, 4, 6)], digits=1))
  range.values <- as.vector(round(summary(colors)[c(1, 6)], digits=1)) - c(0.1, -0.1)
  
  shape.plot <- 
    shape.plot +
    geom_polygon(aes (x = X, y = Y, group = ild, 
                      color = color, fill = color), Q.line.df, size = thickness) +
  geom_text (aes (x = X, y = Y, label = nominhos), color = "black", size = 4) 
  spec.pal <- colorRampPalette (palette, space = 'Lab')
  shape.plot +
    scale_color_gradientn('', colours = spec.pal(10),  breaks = breaks.values, limits =  range.values) +
    scale_fill_gradientn('', colours = spec.pal(10),  breaks = breaks.values, limits =  range.values)
}
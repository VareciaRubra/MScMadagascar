cov.mx<- sp.master.main.data %>% llply(function(x) x$matrix$cov) %>% {. %*% t(.)}
mx.rep<- sp.master.main.data %>% ldply(function(x) x$Mx.Rep$BootsRep[1]) %>% {. %*% t(.)}
n.size<- sp.master.main.data %>% ldply(function(x) x$sample.size) %>% {. %*% t(.)}


mask<- llply(sp.master.main.data,  function(x) !is.na(x$matrix$cov[1]))

stats <- llply(sp.master.main.data, function(x) if (!is.na(x$matrix$cov[1])  ) x$matrix$cov else NA) 
      %>% (., function(x) !is.na(x[1]))
  ldply(MeanMatrixStatistics, .progress = progress_text(char = "."), .inform = T)
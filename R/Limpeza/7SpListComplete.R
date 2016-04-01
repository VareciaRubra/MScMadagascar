
lista.taxo <- sp.main.data %>% 
  llply (function(x) paste0(paste(x$info$Museu, x$info$Tombo, x$info$Sexo, sep = "-" ), "," ) ) %>%  # colando as informações dos bichos e separando com -, depois metendo uma virgula no final
           llply(function (x) gsub("NA", "IND", x)) %>%  #os que nao tem info de sexo trocando pra IND
           llply(function (x) gsub("_", "-", x) ) # substituindo os _ nos numeros de tombo pq o Latex nao gosta


lista.Taxo <- lista.taxo
names(lista.taxo) <- gsub("_", " ", names(lista.Taxo)) 


for (i in 1: length(lista.taxo)) {
  
  lista.Taxo[[i]] <- cat("\\paragraph{\\textit{",names(lista.taxo)[i],"}}", lista.taxo[[i]], "\n" )
}




#Preparando dados para o pedido Gabriel Bridging the gap:
current.data <- sp.main.data
current.data <- gen.main.data

#pegando só os viventes:
mask.extant <- ldply(current.data, function(x) unique(x$info$Status) == "Extant")[,2]
cov.mx <- current.data[mask.extant] %>% llply(function(x) x$matrix$cov)
mask.is.na.cov <- ldply(cov.mx, function(x) is.na(x[1]))[,2]
mask.no.na.cov <- ldply(cov.mx, function(x) !is.na(x[1]))[,2]


#Passo 1: verificar quais sao as especies que tem tamanho amostral muito pequeno para ter matriz
trimed.sp.main.data <- sp.main.data
trimed.sp.main.data

extant.data<- current.data[mask.extant]

#Separar:
#média geométrica por espécie:
gm.mean <- extant.data %>% ldply(function(x) x$gm.mean) 
#médias dos caracteres por especie
ed.means <- extant.data%>% llply(function(x) x$ed.means) 
#Tamanho amostral:
n.size <- extant.data %>% ldply(function(x) x$sample.size) 
#Para os que nao tem matriz, assumir que tamanho amostral é de ZERO e rodar a PhyloW

# PC1 ao PC4 por espécie
PCs1<- current.data %>% ldply(function(x) as.data.frame(eigen(x$matrix$cov)$vectors[,1]) )
PCs2<- current.data %>% ldply(function(x) as.data.frame(eigen(x$matrix$cov)$vectors[,2]) )
PCs3<- current.data %>% ldply(function(x) as.data.frame(eigen(x$matrix$cov)$vectors[,3]) )
PCs4<- current.data %>% ldply(function(x) as.data.frame(eigen(x$matrix$cov)$vectors[,4]) )
#Tamanho amostral:
n.size <- current.data %>% ldply(function(x) x$sample.size) 
#Para os que nao tem matriz, assumir que tamanho amostral é de ZERO e rodar a PhyloW






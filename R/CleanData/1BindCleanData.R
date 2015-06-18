########Colando todos os Clean Data ############
eachPath<-"/home/annapenna/MScMadagascar/Data/"
list.files(path=eachPath, pattern = "([A-z]_Clean+).csv", recursive = T)
filez<-dir(path=eachPath, pattern="([A-z]_Clean+).csv", recursive=TRUE)
#colando o endereço do diretorio inteiro
filez<-paste(eachPath,filez, sep="")
filez
#testando a formatacao de cada arquivo
read.csv(file= filez [1], header=TRUE, sep = ",")
#######################################################################################


###################colando os arquivos de cada colecao ################################
#Lendo todos os arquivos e salvando dentro de um objeto como listas
x<-llply(filez, function(x) read.csv(file= x, header=TRUE, sep = ","))
#install.packages("data.table")
#library(data.table)
#install_github("Rdatatable/data.table", build_vignettes = FALSE)
#library("data.table")
#colando as linhas de acordo com colunas correspondentes
#usando a funcao rbind list do data.table que corresponde a um do.call(rbind...)
#colunas que nao sejam encontradas em alguns arquivos vao retornar NA para os demais "fill = T"
all.brutos.cleaned<- rbindlist(l = x, fill = T)
str(all.brutos.cleaned)
write.csv(all.brutos.cleaned, file = "Data/1Master_Data_Cleaned.csv", row.names=F)

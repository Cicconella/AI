source("https://bioconductor.org/biocLite.R")
#install.packages("oro.dicom")
#install.packages("imager")
#install.packages("mmand")
#install.packages("e1071")
#install.packages("wavelets")
#install.packages("radiomics")
#biocLite("EBImage")

library(oro.dicom) # Leitura do arquivo DICOM
library(imager) # Blur da imagem
library(mmand) # Operacoes morfologicas (abertura e fechamento) 
library(e1071) #Clusterizacao FCM
library(wavelets) # Realizar transformacao DWT
library(radiomics) # Gerar a matrix GLCM e extrair features
library(EBImage) # Identificar Maior Componente de Imagem Binaria


# Plotar Imagens
plot_image <- function(m,nome){
  max_escala = 1#min(max(m)/255,1)
  image(t(m), col=grey(0:64*(max_escala)/64), axes=FALSE, xlab=nome, ylab="")
  png(paste(home,"/plots/",nome,".png",sep=""),width = 512,height = 512)
  image(t(m), col=grey(0:64*(max_escala)/64), axes=FALSE, xlab=nome, ylab="")
  dev.off()
}

plot_matrix <- function(m, nome){
  dx = dim(m)[1]
  dy = dim(m)[2]
  plot_image(m[nrow(m):1,],nome)  
}

limpa_preto <- function(m){
  #dim(m)
  if (length(which((apply(m, 1, sum))==0)) > 0){
    m = m[-which(apply(m, 1, sum)==0),]
  }
  if (length(which(apply(m, 2, sum)==0)) > 0){
    m = m[,-which(apply(m, 2, sum)==0)]
  }
  return(m)
}

# Preprocessamento
normaliza = function(x, lo, hi){
  if(x<lo){
    return(0)
  } else if(x>hi){
    return(255)
  } else{
    x = (255*(x-lo))/(hi-lo)
    return(x)
  }
}

# Preprocessamento
moda <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Limpeza
maior_componente <- function(a) {
  m = bwlabel(a)  
  tam = max(table(m)[-1])
  n = which.max(table(m)[-1])
  m[m!=n]=0
  m[m==n]=1
  return(list(matrix=m, max=tam))
}



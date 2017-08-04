#source("https://bioconductor.org/biocLite.R")
#biocLite("EBImage")

#install.packages("oro.dicom")
#install.packages("oro.nifti")
#install.packages("imager")
#install.packages("mmand")
#install.packages("e1071")
#install.packages("wavelets")
#install.packages("radiomics")
#install.packages("wvtool") #Pacote para a deteccao de bordas canny


library(oro.dicom) # Leitura do arquivo DICOM
library(oro.nifti) # Leitura do arquivo NIFTI
#library(imager) # Blur da imagem
library(mmand) # Operacoes morfologicas (abertura e fechamento) 
library(e1071) #Clusterizacao FCM
library(wavelets) # Realizar transformacao DWT
#library(radiomics) # Gerar a matrix GLCM e extrair features
library(EBImage) # Identificar Maior Componente de Imagem Binaria
#library(wvtool) # Canny deteccao de bordas 

# Plotar Imagens
plota_janelada <- function(m,mini,maxi, nome = "Imagem"){
  EBImage::display( (m-mini)/(maxi-mini),title = nome )
}
plota_imagem <- function(m, nome = "Imagem"){
  mini = min(m)
  maxi = max(m)
  plota_janelada(m,mini,maxi,nome)
  #png(paste(home,"/plots/",nome,".png",sep=""),width = 512,height = 512)
  #image(t(m), col=grey(0:64*(max_escala)/64), axes=FALSE, xlab=nome, ylab="")
  #dev.off()
}
plot_matrix <- function(m, nome){
  plot_image(m[1:nrow(m),],nome)#nrow(m):1,],nome)  
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

## Transformacao Wavelet
dwt_rows <- function(m){
  dim_m = dim(m)[1]
  for (i in 1:dim_m){
    wt = dwt(as.numeric(m[i,]), n.levels=1, filter = "haar")
    m[i,] = c(wt@V$V1,wt@W$W1)
  }
  return(m)
}

dwt_matrix = function(m){
  m = dwt_rows(m)
  m = t(dwt_rows(t(m)))
  return(m)
}

## Transformacao Wavelet Inversa
idwt_rows <- function(m){
  dim_m = dim(m)[1]
  for (i in 1:dim_m){
    linha = m[i,]
    wt = dwt(linha, n.levels=1, filter = "haar")
    wt@V$V1 = as.matrix(linha[1:(dim_m/2)])
    wt@W$W1 = as.matrix(linha[(dim_m/2+1):dim_m])
    m[i,] = idwt(wt)
  }
  return(m)
}

idwt_matrix = function(m){
  m = t(idwt_rows(t(m)))
  m = idwt_rows(m)
  return(m)
}

extrai1Q <- function(m){
  dim_m = dim(m)[1]
  for (i in (dim_m/2+1):dim_m){
    m[i,] = rep(0,dim_m/2)
    m[,i] = rep(0,dim_m/2)
  }
  return(m) 
}

# Encontra posicao x,y a partir dos indices de matrix
desindexa = function(m,n){
  dx = dim(m)[1]
  dy = dim(m)[2]
  x = (n-1)%%dx +1
  y = floor((n-1)/dx)+1
  return(rbind(x,y))
}

# Le o DICOM
le_dicom = function(dir,filename){
  fname <-  paste(dir, filename, sep="/")
  dicom <- readDICOMFile(fname)
  return(dicom$img)
}

##### Normaliza e pre processamento ######
pre_normaliza = function(m){
  #m=ori
  m[m<0] = 0 #remove negativos
  bin = m > otsu(m,range = c(0,max(m)))
  binfill = fillHull(bin)
  binmaior = maior_componente(binfill)$matrix
  fabd = m*binmaior
  #EBImage::display(fabd/max(fabd))
  
  h = density(fabd[fabd>500]) #500 corta os pulmoes
  #plot(h)
  #plot(h, xaxt='n')
  #axis(side=1, at=seq(0,2000, 100), labels=seq(0,2000,100), las=2)
  
  significativo = (h$y[-c(1,2)]>2e-4)
  
  i1 = h$x[which(diff(significativo)==1)]
  i2 = h$x[which(diff(significativo)==-1)]
  lo = min(i1)
  hi = max(i2)
  normalizada = fabd
  for(i in 1:dim(fabd)[1]){
    for(j in 1:dim(fabd)[2]){
      normalizada[i,j] = normaliza(fabd[i,j], lo, hi) 
    }
  }
  return(normalizada)
}


# testes
#m = matrix(c(4,5,2,1,7,4,2,1,rep(5,8),seq(1,8),7,6,6,5,9,8,2,2,
#             4,5,2,1,7,4,2,1,rep(5,8),seq(1,8),7,6,6,5,9,8,2,2),8,byrow=T)
#m2 = dwt_matrix(m)
#m5 = idwt_matrix(m2)
#m3 = extrai1Q(m2)
#m4 = idwt_matrix(m3)

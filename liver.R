#### Libraries and functions ####

#install.packages("oro.dicom")
#biocLite("EBImage")

source("https://bioconductor.org/biocLite.R")
library(oro.dicom) # Leitura do arquivo DICOM
library(imager) # Blur da imagem
library(mmand) # Operacoes morfologicas (abertura e fechamento) 
library(e1071) #Clusterizacao FCM
library(wavelets) # Realizar transformacao DWT
library(radiomics) # Gerar a matrix GLCM e extrair features
library(EBImage) # Identificar Maior Componente de Imagem Binaria

# Plotar Imagens
plot_image <- function(m,nome){
  image(t(m), col=grey(0:64/64), axes=FALSE, xlab=nome, ylab="")
}

plot_matrix <- function(m, nome){
  dx = dim(m)[1]
  dy = dim(m)[2]
  plot_image(m[nrow(m):1,],nome)  
}

# Transformacao Wavelet
dwt_rows <- function(m){
  dim_m = dim(m)[1]
  wt=list()
  for (i in 1:dim_m){
    wt[[i]] = dwt(as.numeric(m[i,]), n.levels=1, filter = "haar")
  }
  m2 = m
  for (i in 1:dim_m){
    v = c(wt[[i]]@V$V1,wt[[i]]@W$W1)
    #print(length(v))
    m2[i,] = v #c(wt[[i]]@V$V1,wt[[i]]@W$W1)
  }
  return(m2)
}

dwt_matrix = function(m){
  m2 = dwt_rows(m)
  m3 = t(dwt_rows(t(m2)))
  return(m3)
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

limpa_preto <- function(m){
  m = m[-which((apply(m, 1, sum))==0),]
  m = m[,-which((apply(m, 2, sum))==0)]
  return(m)
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

##### Nossos dados ###### 

### Para um arquivo ###

#dir = "/media/cicconella/8AA6013CA6012A71/Users/Nina/Documents/Machiron/52490000/52490000/"

#dir = "/Users/ludykong/MaChiron/Data/HCC Lirads 5/"
#dir = "/Users/ludykong/MaChiron/Data/Hemangioma grande lobo esquerdo/"
dir = "/Users/ludykong/MaChiron/Data/52490000/"

#filename = "1.2.840.113619.2.327.3.1091195278.42.1381225064.881.121.dcm"
#filename = "1.2.840.113704.1.111.5852.1422287786.17764.dcm"
filename = 65643294

fname <-  paste(dir, filename, sep="")

abdo <- readDICOMFile(fname)
#names(abdo)

# Dados Importantes do Header

# abdo$hdr
# which(abdo$hdr[,3]=="SliceLocation")
# abdo$hdr[abdo$hdr$name == "StudyDescription",]
# abdo$hdr[abdo$hdr$name == "SeriesDescription",] #30
# abdo$hdr[abdo$hdr$name == "ImageType",] #10 ORIGINAL PRIMARY AXIAL
# abdo$hdr[abdo$hdr$name == "SliceThickness",] #64
# abdo$hdr[abdo$hdr$name == "SliceLocation",] #125
# abdo$hdr[abdo$hdr$name == "PixelSpacing",] #158 0.703125 0.703125


plot_image(abdo$img, "Original")

a = abdo$img
a = as.array(a)
hist(a[a>100], nc=100000)

##### Normalizar #####

lo = 750
hi = 1250

for(i in 1:dim(a)[1]){
  for(j in 1:dim(a)[2]){
   a[i,j] = normaliza(a[i,j], lo, hi) 
  }
}
plot_image(normalizada, "Normalizada")

normalizada = a

##### Recorte da janela com o figado #####  
dim(normalizada)
pikachu = dim(normalizada)[1]
min_x = floor(pikachu*0.75)+1
max_x = pikachu
min_y = 0
max_y = floor(pikachu*0.25)

janela = normalizada
janela = janela[,-c(min_x:max_x)]
janela = janela[-c(min_y:max_y),]
janela = limpa_preto(janela)

plot_image(janela, "Janela")

##### Filtro anisotropico #####

filtrada = as.cimg(janela)
tmp = proc.time() 
filtrada = blur_anisotropic(filtrada,ampl=1e4,sharp=1) 
proc.time()-tmp 

filtrada = as.matrix(filtrada)
plot_image(filtrada, "Filtrada")

##### Analise do histograma #####

hist(filtrada, nc=100)
hist(filtrada[-c(which(filtrada<100), which(filtrada>225))], nc = 1000)
filtrada = round(filtrada)
mod = moda(filtrada[-c(which(filtrada<100), which(filtrada>225))])
mod
dp = 25

##### Aplicacoes morfologicas #####


binaria = filtrada

binaria[binaria>(mod+dp)] = 0
binaria[binaria<(mod-dp)] = 0
binaria[binaria!=0] = 1


plot_image(binaria, "Binaria Pre Morfo")

kernel <- shapeKernel(c(13,13), type="disc")

binaria = opening(binaria, kernel)
binaria = closing(binaria, kernel)

#tamfig = table(binaria)[ ]
plot_image(binaria, "Binaria Pos Morfo")

dim(binaria)
l = maior_componente(binaria)
maior_binaria = l$matrix
tamanho_figado = l$max
plot_image(maior_binaria, "Maior componente da Binaria")

morfo = maior_binaria * filtrada
morfo = limpa_preto(morfo)
plot_image(morfo, "Morfo")

##### FCM #####
# linha = as.vector(morfo)
# #hist(morfo[morfo>100], nc=100)
# 
# resultado = cmeans(linha, iter.max = 150, centers=c(0, mod/2, mod), dist = "manhattan", m = 1.5)
# #(resultado$centers)
# #dim(resultado$membership)
# 
# cluster = matrix(resultado$cluster, nrow = dim(morfo)[1])
# plot_image(cluster,"Clusters")
# 
# kernel = shapeKernel(c(5,5), type="disc")
# 
# limpo = cluster
# limpo = closing(limpo, kernel)
# limpo = opening(limpo, kernel)
# 
# plot_image(limpo, "Limpo Clusters")
# 
# # Selecionar Cluster 2 e 3
# limpo_2 = limpo
# limpo_2[limpo_2 == 1] = 0
# limpo_2[limpo_2 != 0] = 1

##### Extracao de caracteristicas de textura #####
detalhes = morfo

plot_matrix(detalhes, "Lesao")

# Recorta Janela quadrada de lado par
M = min(dim(detalhes)) 
if (M%%2 == 1) M = M-1
detalhes = detalhes[1:M, 1:M]
plot_matrix(detalhes, "Janela Detalhes")

dwt_detalhes=dwt_matrix(detalhes)
plot_matrix(dwt_detalhes, "Detalhes")

# Extracao de features de textura com a matrix de dependencia de niveis de cinza (GLCM)
m = glcm(detalhes, angle=0,d=1)
calc_features(m) #quais features usaremos depende da rede neural

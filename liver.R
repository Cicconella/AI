#### Libraries and functions ####

#install.packages("oro.dicom")
library(oro.dicom) # Leitura do arquivo DICOM
library(imager) # Blur da imagem
library(mmand) # Operacoes morfologicas (abertura e fechamento) 
library(e1071) #Clusterizacao FCM
library(wavelets) # Realizar transformacao DWT
library(radiomics) # Gerar a matrix GLCM e extrair features

plot_image <- function(m,nome){
  image(t(m), col=grey(0:64/64), axes=FALSE, xlab=nome, ylab="")
}

plot_matrix <- function(m, nome){
  dx = dim(m)[1]
  dy = dim(m)[2]
  plot_image(m[nrow(m):1,],nome)  
}

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

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

##### Nossos dados ###### 

### Para um arquivo ###

#dir = "/media/cicconella/8AA6013CA6012A71/Users/Nina/Documents/Machiron/52490000/52490000/"
dir = "/Users/ludykong/MaChiron/Data/HCC Lirads 5/"
filename = "1.2.840.113619.2.327.3.1091195278.42.1381225064.881.121.dcm"
#filename = 65643294
fname <-  paste(dir, filename, sep="")

abdo <- readDICOMFile(fname)
#names(abdo)

#abdo$hdr
abdo$hdr[abdo$hdr$name == "StudyDescription",]
abdo$hdr[abdo$hdr$name == "SeriesDescription",]

#dim(abdo$hdr)
abdo$hdr[30,]
abdo$hdr[12,]
which(abdo$hdr[,3]=="SliceLocation")
abdo$hdr[125,6]


plot_image(abdo$img, "Original")
dim(abdo$img)

a = abdo$img
dim(a)

#min(apply(a, 1, min))
#max(apply(a, 1, max))

a = as.array(a)
length(a)

hist(a[a>100], nc=100000)

sqrt(var(as.vector(a)))

b = as.vector(a[a>750])
b = b[b<1250]

summary(b)
hist(b, nc=100)

##### Normalizar #####

lo = 750
hi = 1250

for(i in 1:dim(a)[1]){
  for(j in 1:dim(a)[2]){
   a[i,j] = normaliza(a[i,j], lo, hi) 
  }
}
plot_image(normalizada, "Normalizada")
#hist(a[a>0])

min(a)
max(a)

dim(a)

normalizada = a
rm(a)

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
janela = janela[,-which((apply(janela, 2, sum))==0)]
janela = janela[-which((apply(janela, 1, sum))==0),]

plot_image(janela, "Janela")

##### Filtro anisotropico #####

filtrada = as.cimg(janela)
tmp = proc.time() 
filtrada = blur_anisotropic(filtrada,ampl=1e4,sharp=1) 
proc.time()-tmp 

filtrada = as.matrix(filtrada)
plot_image(filtrada, "Filtrada")

##### Analise do histograma #####

#hist(filtrada, nc=100)
#hist(filtrada[-c(which(filtrada<50), which(filtrada>225))], nc = 1000)
filtrada = round(filtrada)
moda = Mode(filtrada[-c(which(filtrada<50), which(filtrada>225))])
dp = 25

##### Aplicacoes morofologicas #####

binaria = filtrada

binaria[binaria>(moda+dp)] = 0
binaria[binaria<(moda-dp)] = 0
binaria[binaria!=0] = 1


plot_image(binaria, "Binaria Pre Morfo")

kernel <- shapeKernel(c(13,13), type="disc")

binaria = opening(binaria, kernel)
binaria = closing(binaria, kernel)

plot_image(binaria, "Binaria Pos Morfo")

morfo = binaria * filtrada

#morfo = morfo[-which(apply(morfo, 1, sum)==0),]
#morfo = morfo[,-which(apply(morfo, 2, sum)==0)]
plot_image(morfo, "Morfo")
#hist(morfo[morfo>150], nc = 1000)

##### FCM #####

linha = as.vector(morfo)
length(linha)

#hist(morfo[morfo>100], nc=100)
tmp=proc.time()
resultado = cmeans(linha, iter.max = 150, centers=c(0, moda/2, moda), dist = "manhattan", m = 1.5)
proc.time()-tmp
#(resultado$centers)

#dim(resultado$membership)

#vetor = c(0, 0.85, 0.15)

# clusteriza = function(vetor){
# 
#   if(vetor[1]>0.9){
#     return(1)
#   }else if(vetor[2]>0.5){
#     return(2)
#   } else {
#     return(3)
#   }
#     
# }

#cluster = apply(resultado$membership, 1, clusteriza)

#length(cluster)

#hist(resultado$cluster)

#cluster = matrix(cluster, nrow = 233)

cluster = matrix(resultado$cluster, nrow = dim(morfo)[1])

plot_image(cluster,"Clusters")

kernel = shapeKernel(c(5,5), type="disc")

limpo = cluster

limpo = closing(limpo, kernel)
limpo = opening(limpo, kernel)

plot_image(limpo, "Limpo Clusters")

limpo[limpo == 1] = 0
limpo[limpo != 0] = 1

##### Extracao de caracteristicas de textura #####

detalhes = limpo*morfo

plot_matrix(detalhes, "Lesao")

M = min(dim(detalhes)) 
if (M%%2 == 1) M = M-1
detalhes = detalhes[1:M, 1:M]

plot_matrix(detalhes, "Janela Detalhes")
dwt_detalhes=dwt_matrix(detalhes)
plot_matrix(dwt_detalhes, "Detalhes")

#class(detalhes)
#dim(detalhes)

min(detalhes)
max(detalhes)

# Extracao de features de textura com a matrix de dependencia de niveis de cinza (GLCM)
m = glcm(detalhes, angle=0,d=1)
calc_features(m) #quais features usaremos depende da rede neural
#### Libraries and functions ####

#install.packages("oro.dicom")
library(oro.dicom)
library(imager)
library(mmand)
library(e1071)

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


fname <-  c("/home/cicconella/Desktop/52490000/52490000/65643294")
#fname <-  c("/home/cicconella/Desktop/52490000/52490000/65647276")
#fname <-  c("/home/cicconella/Desktop/52490000/52490000/65646836")
#fname <-  c("/home/cicconella/Desktop/52490000/52490000/65646946")

fname <-  c("/home/cicconella/Desktop/52490000/52490000/65653187")
fname <-  c("/home/cicconella/Desktop/52490000/52490000/65643690")
fname <-  c("/home/cicconella/Desktop/52490000/52490000/65662305")

abdo <- readDICOMFile(fname)
names(abdo)

abdo$hdr

abdo$hdr[abdo$hdr$name == "ContrastBolusStartTime",]
abdo$hdr[abdo$hdr$name == "StudyDescription",]

dim(abdo$hdr)
abdo$hdr[30,]
abdo$hdr[12,]
which(abdo$hdr[,3]=="SliceLocation")
abdo$hdr[110,6]


#png("/home/cicconella/Desktop/52490000/tm1.png")
image(t(abdo$img), col=grey(0:64/64), axes=FALSE, xlab="", ylab="")
#dev.off()

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

#hist(a[a>0])

min(a)
max(a)

dim(a)

normalizada = a

rm(a)

##### Recorte da janela com o figado #####  


janela = normalizada
janela = janela[,-c(256:512)]
janela = janela[-c(0:128),]
janela = janela[,-which((apply(janela, 2, sum))==0)]
janela = janela[-which((apply(janela, 1, sum))==0),]
  
dim(janela)
janela = normalizada
janela = janela[,-c(256:512)]
janela = janela[-c(0:128),]
janela = janela[,-which((apply(janela, 2, sum))==0)]
janela = janela[-which((apply(janela, 1, sum))==0),]
  
dim(janela)
  
    
image(t(janela), col=grey(0:64/64))

##### Filtro anisotropico #####

# png("/home/cicconella/Desktop/52490000/normalizada.png")
#image(t(janela), col=grey(0:64/64), axes=FALSE, xlab="", ylab="")
filtrada = as.cimg(janela)
#plot(filtrada, axes = F, xlab = "", ylab = "")
# dev.off()
tmp = proc.time() 
filtrada = blur_anisotropic(filtrada,ampl=1e4,sharp=1) 
proc.time()-tmp 
#plot(filtrada, axes = F, xlab = "", ylab = "")

filtrada = as.matrix(filtrada)
head(filtrada)
dim(filtrada)

image(t(filtrada), col=grey(0:64/64), axes=FALSE, xlab="", ylab="")

##### Analise do histograma #####

#hist(filtrada, nc=100)
#hist(filtrada[-c(which(filtrada<50), which(filtrada>225))], nc = 1000)

filtrada = round(filtrada)

#image(t(filtrada), col=grey(0:64/64), axes=FALSE, xlab="", ylab="")

moda = Mode(filtrada[-c(which(filtrada<50), which(filtrada>225))])

dp = 25

##### Aplicacoes morofologicas #####

binaria = filtrada

binaria[binaria>(moda+dp)] = 0
binaria[binaria<(moda-dp)] = 0
binaria[binaria!=0] = 1


#image(t(binaria), col=grey(0:64/64))

kernel <- shapeKernel(c(13,13), type="disc")

binaria = opening(binaria, kernel)
binaria = closing(binaria, kernel)

#image(t(binaria), col=grey(0:64/64))
#image(t(filtrada), col=grey(0:64/64))

morfo = binaria * filtrada

#morfo = morfo[-which(apply(morfo, 1, sum)==0),]
#morfo = morfo[,-which(apply(morfo, 2, sum)==0)]

max(morfo)

#image(t(morfo), col=grey(0:(64*194/255)/64))

#hist(morfo, nc = 1000)

#hist(morfo[morfo>150], nc = 1000)

#morfo = morfo[1:121,1:75]

##### FCM #####

linha = as.vector(morfo)
length(linha)

#hist(morfo[morfo>100], nc=100)

# tmp=proc.time()
# resultado = fanny(linha, k=3, metric = "man", cluster.only = TRUE, memb.exp = 1.5)
# proc.time()-tmp

tmp=proc.time()
resultado = cmeans(linha, iter.max = 150, centers=c(0, moda/2, moda), dist = "manhattan", m = 1.5)
proc.time()-tmp

(resultado$centers)

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

image(t(cluster), col=grey(0:64/64))

kernel = shapeKernel(c(5,5), type="disc")

limpo = cluster

limpo = closing(limpo, kernel)
limpo = opening(limpo, kernel)

image(t(limpo), col=grey(0:64/64))



limpo[limpo == 1] = 0
limpo[limpo != 0] = 1



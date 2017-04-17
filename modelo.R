library(imager)
library(e1071)

a = readPNG("Saudavel.png")
a = readPNG("TrianguloMisto.png")
head(a)
class(a)
a = a[,,1]

hist(a)
hist(a[-c(which(a<0.1), which(a>0.3))])
a[a<=0.2] = 0
a[a>0.2 & a < 0.3] = 0.25
a[a>=0.3 & a < 0.6] = 0.4
a[a>=0.6] = 0.8
table(a)
image(a, col=grey(0:64*(max(a))/64), axes=FALSE, ylab="")

linha = as.vector(a)
table(linha)

if(length(table(linha))==2){
  centers = c(0,1)
}else{
  centers = c(0, 0.2, 0.4, 0.8) 
}
  

resultado = cmeans(linha, centers=centers, dist = "manhattan")
(resultado$centers)

cluster = matrix(resultado$cluster, nrow = dim(a)[1])
hist(cluster)
image(cluster, col=grey(0:64/64), axes=FALSE, ylab="")

areas = table(cluster)/sum(table(cluster))
areas


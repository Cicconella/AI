library(imager)
library(e1071)

a = readPNG("Saudavel.png")
a = readPNG("Triangulo.png")
head(a)
class(a)
a = a[,,1]

hist(a)
image(a, col=grey(0:64*(max(a))/64), axes=FALSE, ylab="")

filtrada = as.cimg(a)
filtrada = blur_anisotropic(filtrada,ampl=1e4,sharp=1) 

filtrada = as.matrix(filtrada)
image(filtrada, col=grey(0:64*(max(filtrada))/64), axes=FALSE, ylab="")

linha = as.vector(a)
hist(a, nc=200)

resultado = cmeans(linha, iter.max = 150, centers=3, dist = "manhattan", m = 2)
(resultado$centers)
dim(resultado$membership)

cluster = matrix(resultado$cluster, nrow = dim(a)[1])
image(cluster, axes=FALSE, ylab="")





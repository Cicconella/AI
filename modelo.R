library(png)
library(imager)
library(e1071)

# hist(a)
# hist(a[-c(which(a<0.1), which(a>0.3))])
# a[a<=0.2] = 0
# a[a>0.2 & a < 0.3] = 0.25
# a[a>=0.3 & a < 0.6] = 0.4
# a[a>=0.6] = 0.8

a = readPNG("Saudavel.png")
a = readPNG("Triangulo.png")
a = readPNG("TrianguloBlur.png")
a = a[,,1]

image(a, col=grey(0:64*(max(a))/64), axes=FALSE, ylab="")
linha = as.vector(a)

resultado = cmeans(linha, centers=c(0, 0.25, 0.5, 1)) #, dist = "manhattan")
centers = resultado$centers

cluster = matrix(resultado$cluster, nrow = dim(a)[1])
image(cluster, col=grey(0:64/64), axes=FALSE, ylab="")

areas = table(cluster)/sum(table(cluster))
centers= centers[areas> 0.005]
areas = areas[areas> 0.005]
rbind(centers,areas)

saudavel = which(centers >0.35 & centers<0.45)
if (length(table(areas)) ==2){
  print("esta saudavel")
}else{
  claro = which(centers > 0.45)
  if (length(claro)!=0){
    print(paste0("algo claro de area ",formatC(areas[claro]*100,digits=2),"%"))
  }
  escuro = which(centers > 0.05 & centers < 0.35)
  if (length(escuro)!=0){
    print(paste0("algo escuro de area ",formatC(areas[escuro]*100,digits=2),"%"))
  }
}

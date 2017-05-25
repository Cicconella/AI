install.packages("pnn")
install.packages("neuralnet")

library(png)
library(imager)
library(radiomics)
library(pnn)
library(neuralnet)

normlinha <- function(vetor){
  minimo = min(vetor)
  maximo = max(vetor)
  d = maximo-minimo
  vetor = (vetor - minimo)/d
  return(vetor)
}

normset <- function(dados){
  return(apply(dados, 2, normlinha))
}
  
healthybase = paste(getwd(), "/testimgs/saudavel", sep="")
tribase = paste(getwd(), "/testimgs/triangulo", sep="")
testbase = paste(getwd(), "/testimgs/testes", sep="")
healthyfiles = list.files(healthybase)
trifiles = list.files(tribase)
testfiles = list.files(testbase)

### Gera features saudaveis
featmatrix = c()
for (arq in healthyfiles){
  a = readPNG(paste(healthybase,arq,sep="/"))[,,1]
  #image(a, col=grey(0:64*(max(a))/64), axes=FALSE, ylab="")  
  m = glcm(a, angle=0,d=1)
  f = calc_features(m) #quais features usaremos depende da rede neural
  f=f[1:5]
  featmatrix = rbind(featmatrix,f)
}

### Gera features triangulo
trimatrix = c()
for (arq in trifiles){
  a = readPNG(paste(tribase,arq,sep="/"))[,,1]
  #image(a, col=grey(0:64*(max(a))/64), axes=FALSE, ylab="")  
  m = glcm(a, angle=0,d=1)
  f = calc_features(m) #quais features usaremos depende da rede neural
  f=f[1:5]
  trimatrix = rbind(trimatrix,f)
}
### Gera features test
testmatrix = c()
for (arq in testfiles){
  a = readPNG(paste(testbase,arq,sep="/"))[,,1]
  #image(a, col=grey(0:64*(max(a))/64), axes=FALSE, ylab="")  
  m = glcm(a, angle=0,d=1)
  f = calc_features(m) #quais features usaremos depende da rede neural
  f=f[1:5]
  testmatrix = rbind(testmatrix,f)
}
testmatrix
summary(featmatrix)
summary(trimatrix)
summary(testmatrix)

plot(c(featmatrix[,1],trimatrix[,1]),c(featmatrix[,2],trimatrix[,2]),col = c(rep("blue",6),rep("red",6)),pch = 16)

d = rbind(featmatrix,trimatrix)
d = cbind(c(rep(0,6),rep(1,6)),d)
colnames(d) = c("class","v1","v2","v3","v4","v5")

nd = normset(d)
net.d = neuralnet(class~v1+v2+v3+v4+v5, nd, rep=5, linear.output=FALSE,threshold = 0.001)
plot(net.d,rep="best")
compute(net.d,nd[,2:6])






# Exemplo pnn
# library(pnn)
# data(norms)
# norms
# class(norms)
# apply(d[1,],2,class)
# plot(norms[,2],norms[,3],col = c(rep("blue",200),rep("red",200)),pch = "O")
# pnn <- learn(norms)
# pnn <- smooth(pnn, sigma=0.9)
# pnn = perf(pnn)
# pnn$observed
# pnn$guessed
# pnn$success
# pnn$fails
# pnn$success_rate
# pnn$bic
# ## Not run: pnn <- perf(pnn) # Optional
# ## Not run: pnn$success_rate # Optional
# apply(norms[,2:3], 1, guess,nn=pnn)
# 
# # The short way
# guess(smooth(learn(norms), sigma=0.8), c(1,2.1))
# guess(smooth(learn(norms), sigma=0.8), c(2,1))
# guess(smooth(learn(norms), sigma=0.8), c(1.5,1))
# Demonstrations
## Not run: demo("norms-trainingset", "pnn")
## Not run: demo("small-trainingset", "pnn")

# resultado = cmeans(linha, centers=c(0, 0.25, 0.5, 1)) #, dist = "manhattan")
# centers = resultado$centers
# 
# cluster = matrix(resultado$cluster, nrow = dim(a)[1])
# image(cluster, col=grey(0:64/64), axes=FALSE, ylab="")
# 
# areas = table(cluster)/sum(table(cluster))
# centers= centers[areas> 0.005]
# areas = areas[areas> 0.005]
# rbind(centers,areas)
# 
# saudavel = which(centers >0.35 & centers<0.45)
# if (length(table(areas)) ==2){
#   print("esta saudavel")
# }else{
#   claro = which(centers > 0.45)
#   if (length(claro)!=0){
#     print(paste0("algo claro de area ",formatC(areas[claro]*100,digits=2),"%"))
#   }
#   escuro = which(centers > 0.05 & centers < 0.35)
#   if (length(escuro)!=0){
#     print(paste0("algo escuro de area ",formatC(areas[escuro]*100,digits=2),"%"))
#   }
# }

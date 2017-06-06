#install.packages("pnn")
#install.packages("neuralnet")

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

### Features  
feat = c("glcm_mean","glcm_variance","glcm_energy","glcm_contrast","glcm_entropy","glcm_homogeneity1","glcm_correlation","glcm_IDMN")
length(feat)

healthybase = paste(getwd(), "/testimgs/saudavel", sep="")
tribase = paste(getwd(), "/testimgs/triangulo", sep="")
#testbase = paste(getwd(), "/testimgs/testes", sep="")
healthyfiles = list.files(healthybase)
trifiles = list.files(tribase)
#testfiles = list.files(testbase)

#a = readPNG(paste(healthybase,"s1.png",sep="/"))[,,1]
#image(a, col=grey(0:64*(max(a))/64), axes=FALSE, ylab="")  

### Gera features saudaveis
featmatrix = c()
for (arq in healthyfiles){
  a = readPNG(paste(healthybase,arq,sep="/"))[,,1]
  m = glcm(a, angle=0,d=1)
  f = calc_features(m) #quais features usaremos depende da rede neural
  f = f[names(f)%in% feat]
  featmatrix = rbind(featmatrix,f)
}

### Gera features triangulo
trimatrix = c()
for (arq in trifiles){
  a = readPNG(paste(tribase,arq,sep="/"))[,,1]
  m = glcm(a, angle=0,d=1)
  f = calc_features(m) #quais features usaremos depende da rede neural
  f = f[names(f)%in% feat]
  trimatrix = rbind(trimatrix,f)
}

### Gera features test
# testmatrix = c()
# for (arq in testfiles){
#   a = readPNG(paste(testbase,arq,sep="/"))[,,1]
#   m = glcm(a, angle=0,d=1)
#   f = calc_features(m) #quais features usaremos depende da rede neural
#   testmatrix = rbind(testmatrix,f)
# }

summary(featmatrix)
summary(trimatrix)
apply(featmatrix, 2, mean)
apply(trimatrix, 2, mean)
#summary(testmatrix)

#plot(c(featmatrix[,1],trimatrix[,1]),c(featmatrix[,2],trimatrix[,2]),col = c(rep("blue",6),rep("red",6)),pch = 16)
h_n = dim(featmatrix)[1]
tri_n = dim(trimatrix)[1]
d = rbind(featmatrix,trimatrix)
nd = normset(d)
nd_h = nd[1:h_n,]
nd_t = nd[ (h_n+1):(h_n+tri_n),]

trainp = 0.75
h_train = floor(h_n*trainp)
train_h_set = sample(1:h_n, h_train)
test_h_set = (1:h_n)[-train_h_set]
healthytrain = nd_h[train_h_set,]
healthytest =  nd_h[test_h_set,]
  
tri_n = dim(trimatrix)[1]
tri_train = floor(tri_n*trainp)
train_t_set = sample(1:tri_n, tri_train)
test_t_set = (1:tri_n)[-train_t_set]
tritrain = nd_t[train_t_set,]
tritest =  nd_t[test_t_set,]

train = cbind(c(rep(0,h_train),rep(1,tri_train)),rbind(healthytrain,tritrain) )
colnames(train)[1] = c("class")
x = paste(colnames(train)[-1],collapse="+")
net.d = neuralnet(paste('class ~ ' ,x),train, rep=5, linear.output=FALSE,threshold = 0.001)
plot(net.d,rep="best")
compute(net.d,train[,-1])

test = rbind(healthytest,tritest)
compute(net.d,test)




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

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

### Features para a rede
feat = c("glcm_mean","glcm_variance","glcm_energy","glcm_contrast","glcm_entropy","glcm_homogeneity1","glcm_correlation","glcm_IDMN")
length(feat)

### Leitura dos arquivos
healthybase = paste(getwd(), "/testimgs/saudavel", sep="")
tribase = paste(getwd(), "/testimgs/triangulo", sep="")
healthyfiles = list.files(healthybase)
trifiles = list.files(tribase)

# Para ver os PNGs
#a = readPNG(paste(healthybase,"s1.png",sep="/"))[,,1]
#image(a, col=grey(0:64*(max(a))/64), axes=FALSE, ylab="")  

### Gera as matrizes de features saudaveis e triangulo
featmatrix = c()
for (arq in healthyfiles){
  a = readPNG(paste(healthybase,arq,sep="/"))[,,1]
  m = glcm(a, angle=0,d=1)
  f = calc_features(m) #quais features usaremos depende da rede neural
  f = f[names(f)%in% feat]
  featmatrix = rbind(featmatrix,f)
}

trimatrix = c()
for (arq in trifiles){
  a = readPNG(paste(tribase,arq,sep="/"))[,,1]
  m = glcm(a, angle=0,d=1)
  f = calc_features(m) #quais features usaremos depende da rede neural
  f = f[names(f)%in% feat]
  trimatrix = rbind(trimatrix,f)
}

apply(featmatrix, 2, mean)
apply(trimatrix, 2, mean)

#plot(c(featmatrix[,1],trimatrix[,1]),c(featmatrix[,2],trimatrix[,2]),col = c(rep("blue",6),rep("red",6)),pch = 16)

#### Separa treino e test saudavel e tri
h_n = dim(featmatrix)[1]
tri_n = dim(trimatrix)[1]
trainp = 0.75
h_train = floor(h_n*trainp)
train_h_set = sample(1:h_n, h_train)
test_h_set = (1:h_n)[-train_h_set]
healthytrain = featmatrix[train_h_set,]
healthytest =  featmatrix[test_h_set,]
  
tri_n = dim(trimatrix)[1]
tri_train = floor(tri_n*trainp)
train_t_set = sample(1:tri_n, tri_train)
test_t_set = (1:tri_n)[-train_t_set]
tritrain = trimatrix[train_t_set,]
tritest =  trimatrix[test_t_set,]

### Normaliza todas matrizes
treino = rbind(healthytrain,tritrain)
medias = apply(treino, 2, mean)
desvios = apply(treino, 2, sd)
for( i in 1:length(feat) ){
  treino[,i] = (treino[,i] - medias[i])/desvios[i]
}
teste=rbind(healthytest,tritest)
for( i in 1:length(feat)){
  teste[,i] = (teste[,i] - medias[i])/desvios[i]
}
#d = rbind(featmatrix,trimatrix)
#nd = normset(d)
#nd_h = nd[1:h_n,]
#nd_t = nd[ (h_n+1):(h_n+tri_n),]



##### Treina a rede #####
train = cbind(c(rep(0,h_train),rep(1,tri_train)), treino)
colnames(train)[1] = c("class")
x = paste(colnames(train)[-1],collapse="+")
net.d = neuralnet(paste('class ~ ' ,x),train, rep=5, linear.output=FALSE,threshold = 0.01,act.fct="tanh")
plot(net.d,rep="best")
compute(net.d,train[,-1])

##### Teste #####
compute(net.d,teste)

#### FIM ####
install.packages("pnn")
install.packages("neuralnet")

library(png)
library(imager)
library(radiomics)
library(pnn)
library(neuralnet)

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
# +v2+v3+v4+v5
net.d = neuralnet(class~v1, d, rep=10, linear.output=FALSE)
plot(net.d,rep="best")
compute(net.d,d[,2:2])
############# Neural Net Example #####################
?neuralnet
AND <- c(rep(0,7),1)
OR <- c(0,rep(1,7))
binary.data <- data.frame(expand.grid(c(0,1), c(0,1), c(0,1)), AND, OR)
print(net <- neuralnet(AND+OR~Var1+Var2+Var3, binary.data, hidden=1,
                       rep=10, err.fct="ce", linear.output=FALSE))
plot(net,rep="best")

XOR <- c(0,1,1,0)
xor.data <- data.frame(expand.grid(c(0,1), c(0,1)), XOR)
print(net.xor <- neuralnet(XOR~Var1+Var2, xor.data, hidden=2, rep=5))
plot(net.xor, rep="best")

data(infert, package="datasets")
print(net.infert <- neuralnet(case~parity+induced+spontaneous, infert,
                              err.fct="ce", linear.output=FALSE, likelihood=TRUE))
infertm = t(matrix(c(1,0,0,1,0,1,1,1,1),nrow=3))
compute(net.infert, infertm)
plot(net.infert)
gwplot(net.infert, selected.covariate="parity")
gwplot(net.infert, selected.covariate="induced")
gwplot(net.infert, selected.covariate="spontaneous")
confidence.interval(net.infert)
Var1 <- runif(50, 0, 100)
sqrt.data <- data.frame(Var1, Sqrt=sqrt(Var1))
print(net.sqrt <- neuralnet(Sqrt~Var1, sqrt.data, hidden=10,
                            threshold=0.01))
plot(net.sqrt)
compute(net.sqrt, (1:10)^2)$net.result
Var1 <- rpois(100,0.5)
Var2 <- rbinom(100,2,0.6)
Var3 <- rbinom(100,1,0.5)
SUM <- as.integer(abs(Var1+Var2+Var3+(rnorm(100))))
sum.data <- data.frame(Var1+Var2+Var3, SUM)
print(net.sum <- neuralnet(SUM~Var1+Var2+Var3, sum.data, hidden=1,
                           act.fct="tanh"))
plot(net.sum)

Var1 <- rpois(100,0.5)
Var2 <- rbinom(100,2,0.6)
Var3 <- rbinom(100,1,0.5)
SUM <- as.integer(abs(Var1+Var2+Var3+(rnorm(100))))
sum.data <- data.frame(Var1+Var2+Var3, SUM)
print(net.sum <- neuralnet( SUM~Var1+Var2+Var3,  sum.data, hidden=1, 
                            act.fct="tanh"))
main <- glm(SUM~Var1+Var2+Var3, sum.data, family=poisson())
full <- glm(SUM~Var1*Var2*Var3, sum.data, family=poisson())
prediction(net.sum, list.glm=list(main=main, full=full))













library(datasets)
data("infert")
infert
head(infert)
inf2=infert[,-1]
head(inf2)
summary(infert)
pnn = learn(inf2,category.column = 1)
pnn <- smooth(pnn, sigma=0.9)
pnn <- perf(pnn)
guess(pnn, c(76.04,6,1,2,1,3) )

pnn = learn(d)
pnn <- smooth(pnn, sigma=0.9)
guess(nn=pnn, c(10.926001,209.26467,326.93921,978949.1,13509.294) )
apply(d2, 1, guess,nn=pnn)

guess(nn=pnn, c(10.241536,61.90144,98.63208,108773.8,1952.640) )
guess(smooth(learn(d), sigma=0.8), c(d[1,2:6]))
d2 = rbind(d[,2:6],testmatrix)
d2


guess(pnn,c(0,2,3,4,5))

class()
class(testmatrix[])
guess(pnn,testmatrix[1,])








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

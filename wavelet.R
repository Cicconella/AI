#install.packages("wavelets")
require("wavelets")
library("jpeg")

plot_matrix <- function(m){
  dx = dim(m)[1]
  dy = dim(m)[2]
  image(t(m[nrow(m):1,]),col=grey(0:64/64), axes=FALSE,xlab="", ylab="")  
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

img = readJPEG("house.jpeg")
m1 = img[,,1]
#dim(m2)
m2= m1[41:168,41:168]

#image(m2,col=grey(0:64/64), axes=FALSE, xlab="", ylab="")
plot_matrix(m2)
m2t=dwt_matrix(m2)
plot_matrix(m2t)

# Teste com numeros
m = t(matrix(1:16,nrow=4))
mt=dwt_matrix(m)

print("Oi N")
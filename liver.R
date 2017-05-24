#### Libraries and functions ####

#install.packages("oro.dicom")
#biocLite("EBImage")

source("library.R")

home = getwd()# "/Users/ludykong/MaChiron/MaChironGit"

# Transformacao Wavelet
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

# Preprocessamento
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

# Preprocessamento
moda <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Limpeza
maior_componente <- function(a) {
  m = bwlabel(a)  
  tam = max(table(m)[-1])
  n = which.max(table(m)[-1])
  m[m!=n]=0
  m[m==n]=1
  return(list(matrix=m, max=tam))
}

#### Calcula matriz onde cada elemento e a media dos elementos vizinhos na original
matriz_media = function(m){
  mr = m
  dx = dim(m)[1]
  dy = dim(m)[2]
  for(i in 1:dx){
    for (j in 1:dy){
      v = c(m[i,j])
      if (i>1) v= c(v,m[i-1,j])
      if (i<dx) v= c(v,m[i+1,j])
      if (j>1) v= c(v,m[i,j-1])
      if (j<dy) v= c(v,m[i,j+1])
      mr[i,j] = mean(v)
    }
  }
  return (mr)
}


##### Nossos dados ###### 

### Para um arquivo ###

#dir = "/media/cicconella/8AA6013CA6012A71/Users/Nina/Documents/Machiron/52490000/52490000/"

#dir = "/Users/ludykong/MaChiron/Data/HCC Lirads 5/"
#dir = "/Users/ludykong/MaChiron/Data/Hemangioma grande lobo esquerdo/"
dir = "/Users/ludykong/MaChiron/Data/52490000/"

dir = "/home/cicconella/"

# a = format(1, width = 2, zero.print = T)

#filename = "1.2.840.113619.2.327.3.1091195278.42.1381225064.881.121.dcm"
#filename = "1.2.840.113704.1.111.5852.1422287786.17764.dcm"
filename = 65643294

fname <-  paste(dir, filename, sep="")

abdo <- readDICOMFile(fname)
dim(abdo$hdr)
print(abdo$hdr)
write.table(abdo$hdr,"header")
plot_image(abdo$img, "Original")

a = abdo$img
a = as.array(a)
hist(a[a>10], nc=1000,main = "Histograma Original")

##### Normalizar #####
lo = 800
hi = 1200

normalizada = a
for(i in 1:dim(a)[1]){
  for(j in 1:dim(a)[2]){
    normalizada[i,j] = normaliza(a[i,j], lo, hi) 
  }
}

plot_image(normalizada, "Normalizada")
hist(normalizada[normalizada>10], nc=1000,main = "Histograma Normalizada")

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
janela = limpa_preto(janela)

plot_image(janela, "Janela")

##### Filtro anisotropico #####
filtrada = as.cimg(janela)
tmp = proc.time() 
filtrada = blur_anisotropic(filtrada,ampl=1e4,sharp=1) 
proc.time()-tmp 

filtrada = as.matrix(filtrada)
plot_image(filtrada, "Filtrada")

##### Analise do histograma #####
lim_inf = 50
lim_sup = 200
#dp = 35
length(unique(as.vector(filtrada)))
filtrada = round(filtrada)
hist(filtrada[-c(which(filtrada<lim_inf), which(filtrada>lim_sup))],nc = 256)#[-c(which(filtrada<lim_inf), which(filtrada>lim_sup))], nc = 256)
h = density(filtrada[-c(which(filtrada<lim_inf), which(filtrada>lim_sup))])
plot(h)
h$x
h$y

picos <- h$x[which(diff(sign(diff(h$y )))==-2)]
vales <- h$x[which(diff(sign(diff(h$y )))==2)]

np = length(picos)
np
V1= vales[np-1]
V2= picos[np]

limite_y = h$y[which(h$x == V2)]/30
V3 = h$x[-c(which(h$y > limite_y), which(h$x< V2) )][1]
c(V1,V2,V3)
#mod = moda(filtrada[-c(which(filtrada<lim_inf), which(filtrada>lim_sup))])

##### Aplicacoes morfologicas #####
binaria = filtrada

binaria[binaria>V3] = 0
binaria[binaria<V1] = 0
binaria[binaria!=0] = 1

plot_image(binaria, "Binaria Pre Morfo")

kernel <- shapeKernel(c(5,5), type="disc")
binaria_pos = opening(binaria, kernel)
binaria_pos = closing(binaria_pos, kernel)

#binaria_pos = binaria
plot_image(binaria_pos, "Binaria Pos Morfo")

###### Encontrar maior Componente da Binaria ######
dim(binaria_pos)
l = maior_componente(binaria_pos)
maior_binaria = l$matrix
tamanho_figado = l$max
plot_image(maior_binaria, "Maior componente da Binaria")

masc = fillHull(maior_binaria)
plot_image(masc, "Masc")

morfo = masc * filtrada
#morfo = limpa_preto(morfo)
plot_image(morfo, "Morfo pre morfo")

kernel <- shapeKernel(c(7,7), type="disc")
masc = opening(masc, kernel)
masc = closing(masc, kernel)
plot_image(masc, "Masc pos morfo")

morfo = masc * filtrada
#morfo = limpa_preto(morfo)
plot_image(morfo, "Morfo")

##### FCM #####
#media_morfo = matriz_media(morfo)
#plot_image(media_morfo, "Media Morfo")

#linha = as.vector(media_morfo)
#hist(morfo[morfo>100], nc=200)

#resultado = cmeans(linha, iter.max = 150, centers=c(0, mod-dp, mod), dist = "manhattan", m = 1.5)
#(resultado$centers)
#dim(resultado$membership)

#cluster = matrix(resultado$cluster, nrow = dim(morfo)[1])
#plot_image(cluster,"Clusters")
# 
# kernel = shapeKernel(c(5,5), type="disc")
# 
# limpo = cluster
# limpo = closing(limpo, kernel)
# limpo = opening(limpo, kernel)
# 
# plot_image(limpo, "Limpo Clusters")
# 
# # Selecionar Cluster 2 e 3
# limpo_2 = limpo
# limpo_2[limpo_2 == 1] = 0
# limpo_2[limpo_2 != 0] = 1

##### Extracao de caracteristicas de textura #####
detalhes = morfo

plot_matrix(detalhes, "Lesao")

# Recorta Janela quadrada de lado par
M = min(dim(detalhes)) 
if (M%%2 == 1) M = M-1
detalhes = detalhes[1:M, 1:M]
plot_matrix(detalhes, "Janela Detalhes")

dwt_detalhes=dwt_matrix(detalhes)
plot_matrix(dwt_detalhes, "Detalhes")

# Extracao de features de textura com a matrix de dependencia de niveis de cinza (GLCM)
m = glcm(detalhes, angle=0,d=1)
calc_features(m) #quais features usaremos depende da rede neural


TESTE
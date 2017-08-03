#### Libraries and functions ####
source("library.R")

#home = getwd()# "/Users/ludykong/MaChiron/MaChironGit"

#### Calcula matriz onde cada elemento e a media dos elementos vizinhos na original
# matriz_media = function(m){
#   mr = m
#   dx = dim(m)[1]
#   dy = dim(m)[2]
#   for(i in 1:dx){
#     for (j in 1:dy){
#       v = c(m[i,j])
#       if (i>1) v= c(v,m[i-1,j])
#       if (i<dx) v= c(v,m[i+1,j])
#       if (j>1) v= c(v,m[i,j-1])
#       if (j<dy) v= c(v,m[i,j+1])
#       mr[i,j] = mean(v)
#     }
#   }
#   return (mr)
# }

##### Nossos dados ###### 
#dir da ana
dir = "/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/"

#dir do luis
#dir = "/Users/ludykong/GDrive/"

filename = "MaChiron/Exames/1775933/DICOM/ARTERIAL/12865.dcm"
#files = list.files(dir)
#if (files[length(files)] == "Icon\r"){
#  files = files[-length(files)]
#}
#dir = "/Users/ludykong/MaChiron/Data/52490000/"

fname <-  paste(dir, filename, sep="")

abdo <- readDICOMFile(fname)
ori = abdo$img
#plot_image(ori, "Original")
plota_imagem(ori)
#hist(ori[ori>10], nc=1000,main = "Histograma Original")
##### Normalizar #####
#lo = 800
#hi = 1200
#normalizada = a
#for(i in 1:nrow(a)){
#  for(j in 1:ncol(a)){
#    normalizada[i,j] = normaliza(a[i,j], lo, hi) 
#  }
#}
norm = pre_normaliza(ori)
plota_imagem(norm)

#hist(normalizada[normalizada>10], nc=1000,main = "Histograma Normalizada")

##### Recorte da janela com o figado #####  
min_x = floor(nrow(norm)*0.7)+1
max_x = nrow(norm)
min_y = 0
max_y = floor(nrow(norm)*0.25)

janela = norm
janela = janela[,-c(min_x:max_x)]
janela = janela[-c(min_y:max_y),]
janela = limpa_preto(janela)

plota_imagem(janela)

##### Filtro anisotropico #####
#filtrada = as.cimg(janela)
#tmp = proc.time() 
#filtrada = blur_anisotropic(filtrada,ampl=1e4,sharp=1) 
#proc.time()-tmp 
#filtrada = as.matrix(filtrada)
#plot_image(filtrada, "Filtrada")
#EBImage::display(filtrada/max(filtrada))

##### Analise do histograma #####
lim_inf = 50
lim_sup = 200
#dp = 35
#norm = round(norm)
hist(janela[-c(which(janela<lim_inf), which(janela>lim_sup))],nc = 256)#[-c(which(filtrada<lim_inf), which(filtrada>lim_sup))], nc = 256)
h = density(janela[-c(which(janela<lim_inf), which(janela>lim_sup))])
plot(h)
picos <- h$x[which(diff(sign(diff(h$y )))==-2)]
vales <- h$x[which(diff(sign(diff(h$y )))==2)]

np = length(picos)
np
V1= vales[np-1]
V2= picos[np]

limite_y = h$y[which(h$x == V2)]/30
V3 = h$x[-c(which(h$y > limite_y), which(h$x< V2) )][1]


limite_y = h$y[which(h$x == V2)]/2
mX = h$x[-c(which(h$y > limite_y), which(h$x> V2) )]
V4 = mX[length(mX)]
c(V1,V2,V3,V4)
V1 = V4
#mod = moda(filtrada[-c(which(filtrada<lim_inf), which(filtrada>lim_sup))])

##### Aplicacoes morfologicas #####
binaria = janela

binaria[binaria>V3] = 0
binaria[binaria<V1] = 0
binaria[binaria!=0] = 1

plota_imagem(binaria)

kernel <- shapeKernel(c(3,3), type="disc")
binaria_pos = opening(binaria, kernel)
binaria_pos = closing(binaria_pos, kernel)

#binaria_pos = binaria
#plot_image(binaria_pos, "Binaria Pos Morfo")
plota_imagem(binaria_pos)

###### Encontrar maior Componente da Binaria ######
l = maior_componente(binaria_pos)
maior_binaria = l$matrix
tamanho_figado = l$max
#plot_image(maior_binaria, "Maior componente da Binaria")
plota_imagem(maior_binaria)

masc = fillHull(maior_binaria)
plota_imagem(masc)
#plot_image(masc, "Figado Bin I")

morfo = masc * janela
#morfo = limpa_preto(morfo)
plota_imagem(morfo)
#plot_image(morfo, "Figado I")
kernel <- shapeKernel(c(7,7), type="disc")
masc = closing(masc, kernel)
masc = opening(masc, kernel)
plota_imagem(masc)
morfo = masc * janela
plota_imagem(morfo/max(morfo))
#morfo = limpa_preto(morfo)
#plot_image(morfo, "Morfo")

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

#plot_matrix(detalhes, "Lesao")
plota_imagem(detalhes/max(detalhes))

# Recorta Janela quadrada de lado par
M = min(dim(detalhes)) 
if (M%%2 == 1) M = M-1
detalhes = detalhes[1:M, 1:M]
#plot_matrix(detalhes, "Janela Detalhes")

dwt_detalhes=dwt_matrix(detalhes)
dwt1 = dwt_detalhes[1:M/2,1:M/2]
EBImage::display(dwt1/max(dwt1))
#plot_matrix(dwt_detalhes, "Detalhes")
#plot_matrix(dwt1, "Detalhes 1Q")
# Extracao de features de textura com a matrix de dependencia de niveis de cinza (GLCM)
#m = glcm(detalhes, angle=0,d=1)
#calc_features(m) #quais features usaremos depende da rede neural


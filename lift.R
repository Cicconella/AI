##### Bibliotecas #####
source("library.R")

##### Ler Nifti###### 
#dir da ana
#dir = "/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/"

#dir do luis
dir = "/Users/ludykong/GDrive/"

N_nii = 6
N_img = 460
filename = paste0("MaChiron/Exames/TreinoLITS/volume-",N_nii,".nii")

fname <-  paste(dir, filename, sep="")

# So uma vez
abdo <- readNIfTI(fname)
pixel_data = pixdim(abdo)
abdo_data = img_data(abdo)
dim(abdo_data)
rm(abdo)
gc(reset=T)

##### Seleciona uma faixa #####
<<<<<<< HEAD
m = abdo_data[,,N_img]
=======
m = img_data(abdo)

m = m[,,50]
plota_imagem(m)

m = m[,,N_img]
>>>>>>> c371166e3361cb5ece81155232de7a203cd22f3c
plota_imagem(m,"Original")

m = m + abs(min(m[m!=min(m)]))
#plota_imagem(m,"Sem Negativos")

norm = pre_normaliza(m)

plota_imagem(norm)

plota_imagem(norm,"Normalizada")



##### Recorte da janela com o figado #####  
janela = norm
#plota_imagem(janela, "Janela")

##### Analise do histograma #####
lim_inf = 50
lim_sup = 200
#hist(janela[-c(which(janela<lim_inf), which(janela>lim_sup))], col="blue")#[-c(which(filtrada<lim_inf), which(filtrada>lim_sup))], nc = 256)
h = density(janela[-c(which(janela<lim_inf), which(janela>lim_sup))])
#plot(h)
picos <- h$x[which(diff(sign(diff(h$y )))==-2)]
vales <- h$x[which(diff(sign(diff(h$y )))==2)]

np = length(picos)
V1= vales[np-1]
V2= picos[np]

limite_y = h$y[which(h$x == V2)]/30
V3 = h$x[-c(which(h$y > limite_y), which(h$x< V2) )][1]

limite_y = h$y[which(h$x == V2)]/2
mX = h$x[-c(which(h$y > limite_y), which(h$x> V2) )]
V4 = mX[length(mX)]
#c(V1,V2,V3,V4)

##### Aplicacoes morfologicas #####
binaria = janela

binaria[binaria>V3] = 0
binaria[binaria<V1] = 0
binaria[binaria!=0] = 1

#plota_imagem(binaria, "1a binaria")
kernel <- shapeKernel(c(3,3), type="disc")
binaria_pos = closing(binaria, kernel)
binaria_pos = opening(binaria_pos, kernel)
#plota_imagem(binaria_pos, "1a binaria+OM")

###### Encontrar maior Componente da Binaria ######
l = maior_componente(binaria_pos)
maior_binaria = l$matrix
tamanho_figado = l$max
#plota_imagem(maior_binaria, "Maior_Comp")

masc = fillHull(maior_binaria)
<<<<<<< HEAD
plota_imagem(masc)

morfo = masc * janela
#plota_imagem(morfo)
#plota_imagem(masc, "maior_Comp_Filled")

morfo = masc * janela
#plota_imagem(morfo,"Morfo sem OM final")


m=masc
masc=m
kernel <- shapeKernel(c(7,7), type="disc")
masc = closing(masc, kernel)
masc = opening(masc, kernel)
#plota_imagem(masc,"Mascara pos OM")

morfo = masc * janela

#plota_imagem(morfo)

#plota_imagem(maior_binaria)
#plota_imagem(masc)

##### Comparar o gabarito LiTS #####
filename = "MaChiron/Exames/Teste LITS/segmentation-28.nii"

plota_imagem(morfo,"Morfo pos OM")

##### Comparar o gabarito LiTS #####
<<<<<<< HEAD
filename = paste0("MaChiron/Exames/TreinoLITS/segmentation-",N_nii,".nii")
=======
filename = "MaChiron/Exames/Teste LITS/segmentation-0.nii"

>>>>>>> c371166e3361cb5ece81155232de7a203cd22f3c
fname <-  paste(dir, filename, sep="")

#Rodar so 1 vez por exame
seg <- readNIfTI(fname)
seg_data = img_data(seg)

<<<<<<< HEAD
rm(seg)
gc(reset=T)

seg_m = seg_data[,,N_img]
plota_imagem(seg_m,"Gabarito LiTS Original")

table(seg_m)
table(masc == seg_m)
MASC1 = which(masc==1)
MASC0 = which(masc==0)
GABA1 = which(seg_m==1 | seg_m==2)
GABA0 = which(seg_m==0)

TP = length(intersect(MASC1,GABA1))
TN = length(intersect(MASC0,GABA0))
FP = length(intersect(MASC1,GABA0))
FN = length(intersect(MASC0,GABA1))
PRE = TP/(TP+FP)
SENSI = TP/(TP+FN)
ESPE = TN/(TN+FP)
DICE = 2*TP/(2*TP + FP + FN)
JAC = TP/(TP + FP + FN)
=======
m = img_data(abdo)
m = m[,,60]
plota_imagem(m)
plota_imagem(masc)



seg_m = img_data(seg)
seg_m = seg_m[,,N_img]
plota_imagem(seg_m,"Gabarito LiTS Original")

>>>>>>> c371166e3361cb5ece81155232de7a203cd22f3c

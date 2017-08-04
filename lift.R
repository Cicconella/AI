##### Bibliotecas #####
source("library.R")

##### Ler Nifti###### 
#dir da ana
dir = "/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/"

#dir do luis
#dir = "/Users/ludykong/GDrive/"

filename = "MaChiron/Exames/Teste LITS/volume-0.nii"

fname <-  paste(dir, filename, sep="")

abdo <- readNIfTI(fname)
#abdo
#image(abdo,plot.type="single",z=60)

##### Seleciona uma faixa #####
m = img_data(abdo)
m = m[,,50]
plota_imagem(m)

m = m + abs(min(m[m!=min(m)]))
#plota_imagem(m)

norm = pre_normaliza(m)
plota_imagem(norm)


##### Recorte da janela com o figado #####  
janela = norm

#plota_imagem(janela)

##### Analise do histograma #####
lim_inf = 50
lim_sup = 200
#dp = 35
#norm = round(norm)
hist(janela[-c(which(janela<lim_inf), which(janela>lim_sup))],nc = 256)#[-c(which(filtrada<lim_inf), which(filtrada>lim_sup))], nc = 256)
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

#plota_imagem(binaria)

kernel <- shapeKernel(c(3,3), type="disc")
binaria_pos = opening(binaria, kernel)
binaria_pos = closing(binaria_pos, kernel)

#plota_imagem(binaria_pos)

###### Encontrar maior Componente da Binaria ######
l = maior_componente(binaria_pos)
maior_binaria = l$matrix
tamanho_figado = l$max
#plota_imagem(maior_binaria)

masc = fillHull(maior_binaria)
plota_imagem(masc)

morfo = masc * janela
#plota_imagem(morfo)

m=masc
masc=m
kernel <- shapeKernel(c(7,7), type="disc")
masc = closing(masc, kernel)
masc = opening(masc, kernel)
plota_imagem(masc)

morfo = masc * janela
#plota_imagem(morfo)

#plota_imagem(maior_binaria)
#plota_imagem(masc)

##### Comparar o gabarito LiTS #####
filename = "MaChiron/Exames/Teste LITS/segmentation-28.nii"

fname <-  paste(dir, filename, sep="")

abdo <- readNIfTI(fname)
#image(abdo,plot.type="single",z=60)

m = img_data(abdo)
m = m[,,60]
plota_imagem(m)
plota_imagem(masc)



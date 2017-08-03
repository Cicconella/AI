source("library.R")

home = getwd() # "/Users/ludykong/MaChiron/MaChironGit"
dir = "/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/"
dir = "/Users/ludykong/GDrive/"

#filename = 65643690#65643294

filename = "MaChiron/Exames/1775933/DICOM/ARTERIAL/12760.dcm"
fname <-  paste(dir, filename, sep="")

dicom <- readDICOMFile(fname)
plot_imagem(dicom$img, "Original")
ori = dicom$img
hdr = dicom$hdr
##### Normalizar #####
lo = 800
hi = 1200

normalizada = ori
for(i in 1:dim(ori)[1]){
  for(j in 1:dim(ori)[2]){
    normalizada[i,j] = normaliza(ori[i,j], lo, hi) 
  }
}

plot_image(normalizada, "Normalizada")
hist(normalizada[normalizada>10], nc=1000,main = "Histograma Normalizada")

##### Recorte da janela com o baco #####  
dim(normalizada)
npixels = dim(normalizada)[1]
min_x = floor(npixels*0.6)
max_x = floor(0.8*npixels)
min_y = floor(npixels*0.35)+1
max_y = floor(npixels*0.55)+1

janela = normalizada
janela = janela[c(min_y:max_y),c(min_x:max_x)]
janela = limpa_preto(janela)
janela
plot_image(janela, "Janela")

##### Filtro anisotropico #####
filtrada = as.cimg(janela)
tmp = proc.time() 
filtrada = blur_anisotropic(filtrada,ampl=1e4,sharp=1) 
proc.time()-tmp 

filtrada = as.matrix(filtrada)
plot_image(filtrada, "Filtrada")

##### Analise do histograma #####
lim_inf = 10
lim_sup = 250
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
picos
vales
V1= vales[np-2]
V2= picos[np-1]

limite_y = h$y[which(h$x == V2)]/10
V3 = h$x[-c(which(h$y > limite_y), which(h$x< V2) )][1]
c(V1,V2,V3)
#V1= 125
#V3= 200
#mod = moda(filtrada[-c(which(filtrada<lim_inf), which(filtrada>lim_sup))])

binaria = filtrada
binaria[binaria>V3] = 0
binaria[binaria<V1] = 0
binaria[binaria!=0] = 1

plot_image(binaria, "Binaria Pre Morfo")

##### Aplicacoes morfologicas #####
kernel <- shapeKernel(c(5,5), type="disc")
binaria_pos = opening(binaria, kernel)
binaria_pos = closing(binaria_pos, kernel)
plot_image(binaria_pos, "Binaria Pos Morfo")

###### Encontrar maior Componente da Binaria ######
l = maior_componente(binaria_pos)
maior_binaria = l$matrix
tamanho_baco = l$max
plot_image(maior_binaria, "Maior componente da Binaria")
print(paste("Tamanho do Baco:",tamanho_baco))

masc = fillHull(maior_binaria)
plot_image(masc, "Maior componente Preenchido")

morfo = masc * filtrada
#morfo = limpa_preto(morfo)
plot_image(morfo, "Mascara pre morfo")

kernel <- shapeKernel(c(7,7), type="disc")
masc = opening(masc, kernel)
masc = closing(masc, kernel)
plot_image(masc, "Mascara pos morfo")

morfo = masc * filtrada
#morfo = limpa_preto(morfo)
plot_image(morfo, "Morfo")

a = which(hdr[,3]== "PixelSpacing")
pixelSpacing = hdr[a,6]
b = unlist(strsplit(pixelSpacing, split = " "))
pixel_x = b[1]
pixel_y = b[2]

area_baco = tamanho_baco*as.numeric(pixel_x)*as.numeric(pixel_y)
tamanho_baco
area_baco

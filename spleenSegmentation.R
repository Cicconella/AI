source("library.R")

home = getwd()# "/Users/ludykong/MaChiron/MaChironGit"
dir = "/Users/ludykong/MaChiron/Data/52490000/"
dir = "/Users/ludykong/MaChiron/Data/HCC Lirads 4/"
#list.files(dir)
#filename = 65643294
#Lirads4 
#filename = "1.2.840.113619.2.327.3.1091195278.193.1456136545.506.24.dcm"
#fname <-  paste(dir, filename, sep="")
fname = "/Users/ludykong/GDrive/MaChiron/Exames/Cisto/DICOM/VOLUME ABD SC 1.0  B20f/0006-65662095"
dicom <- readDICOMFile(fname)
plot_image(dicom$img, "Original")
ori = dicom$img
hdr = dicom$hdr

##### Pre Processamento ######
display(ori/max(ori))
bin = ori > otsu(ori,range = c(0,max(ori)))
display(bin)
binfill = fillHull(bin)
display(binfill)
binmaior = maior_componente(binfill)$matrix
display(binmaior)
fabd = ori*binmaior
display(fabd/max(fabd))
plot_image(fabd, "PreProc")
#se quiser pode usar filtro gaussiano y = gblur(x, sigma=1)

##### Normalizar #####
hist(fabd[fabd>1],nc=1000)
h = density(fabd[fabd>100])
plot(h)
#+1 para corrigir o deslocamento do diff
significativo = (h$y[-c(1,2)]>1e-4)
maximoLocal = diff(sign(diff(h$y) ) )==-2
picos = h$x[ (which(significativo & maximoLocal) )+1 ] 
p1 = picos[1]
p2 = picos[2]

#### Encontrar intervalo de valores significativos #
i1 = h$x[which(diff(significativo)==1)]
i2 = h$x[which(diff(significativo)==-1)]

#Hfpre = table(fabd)[-1]
#f1 = as.numeric(names(which(Hfpre == max(Hfpre))))
#Hfpre[names(Hfpre) == f1]
#Hfpre2 = table(fabd[fabd>f1])
#f2 = as.numeric(names(which(Hfpre2 == max(Hfpre2))))
#Hfpre2[names(Hfpre2) == f2]
lo = i1
hi = i2
normalizada = fabd
for(i in 1:dim(fabd)[1]){
  for(j in 1:dim(fabd)[2]){
    normalizada[i,j] = normaliza(fabd[i,j], lo, hi) 
  }
}
display(normalizada/max(normalizada))
plot_image(normalizada, "Normalizada")
hist(normalizada[normalizada>1], nc=1000,main = "Histograma Normalizada")

# Reconstroi usando so a aproximada da DWT
wave = dwt_matrix(normalizada)
display(wave/max(wave))
wave = extrai1Q(wave)
display(wave/max(wave))
wave2 = idwt_matrix(wave)
display(wave2/max(wave2))

# Novo otsu
bin = wave2 > otsu(wave2,range = c(0,max(wave2)))
fB = wave2*bin
display(fB/max(fB))

muw = mean(fB[which(fB != 0)])
sigw = sd(fB[which(fB != 0)])
tw = muw + sigw *2

bin = fB>tw
display(bin)

bin = fillHull(bin)
display(bin)
ff = fB*bin
display(ff/max(ff))


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
?fillHull

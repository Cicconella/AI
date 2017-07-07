source("library.R")

##### Normaliza e pre processamento ######
pre_normaliza = function(m){
  #m=ori
  m[m<0] = 0 #remove negativos
  bin = m > otsu(m,range = c(0,max(m)))
  binfill = fillHull(bin)
  binmaior = maior_componente(binfill)$matrix
  fabd = m*binmaior
  #EBImage::display(fabd/max(fabd))
  
  h = density(fabd[fabd>500]) #500 corta os pulmoes
  #plot(h)
  plot(h, xaxt='n')
  axis(side=1, at=seq(0,2000, 100), labels=seq(0,2000,100), las=2)
  
  significativo = (h$y[-c(1,2)]>2e-4)
  
  i1 = h$x[which(diff(significativo)==1)]
  i2 = h$x[which(diff(significativo)==-1)]
  lo = min(i1)
  hi = max(i2)
  normalizada = fabd
  for(i in 1:dim(fabd)[1]){
    for(j in 1:dim(fabd)[2]){
      normalizada[i,j] = normaliza(fabd[i,j], lo, hi) 
    }
  }
  return(normalizada)
}

#### Encontra ossos e objetos brancos #####
acha_brancos = function(m){
  #bin = m > otsu(m,range = c(0,max(m)))
  #fB = m*bin
  #muw = mean(fB[which(fB != 0)])
  #sigw = sd(fB[which(fB != 0)])
  #tw = muw + sigw
  #bin = fB>tw
  #outro jeito
  bin = matrix(data = 0, ncol = ncol(m), nrow = nrow(m) )
  bin[m==max(m)] = 1
  bin = fillHull(bin)
  #EBImage::display(bin)
  #kernel <- shapeKernel(c(3,3), type="disc")
  #bin_pos = opening(bin, kernel)
  #bin_pos = closing(bin_pos, kernel)
  #EBImage::display(bin_pos)
  # Remover componentes pequenos da mascara
  limiteTamanho = 30
  cc = bwlabel(bin)  
  
  rem = as.numeric(names(which(table(cc) < limiteTamanho)))
  for(i in rem){
    cc[cc == i] = 0
  }
  cc[ cc >1]= 1
#  EBImage::display(cc)
  return (cc)
}

#### Encontra a janela de interesse #####
encontraROI = function(m,vert){
  ## Detectando a ROI
  l = which(apply(m,1,sum) > 0)
  comAbd = l[1]
  fimAbd = l[length(l)]
  l = which(apply(vert,1,sum) > 0)
  comVert = l[1]
  d = comVert - comAbd
  
  l = which(apply(m,2,sum) > 0)
  comAbdv = l[1]
  fimAbdv = l[length(l)]
  meioAbdv = round((comAbdv +fimAbdv)/2)
  #ini_x,fim_x,ini_y fim_y
  return(c(comVert,(fimAbd-d),meioAbdv,fimAbdv))
}

#### Remove componentes escuros e pequenos #####
filtra_componentes = function(ROI, startslice){
  mur = mean(ROI[ROI>0])
  sigr = sd(ROI[ROI>0]) 
  tr = mur + sigr
  
  binROI = ROI
  binROI[binROI < tr] = 0
  binROI[binROI >= tr] = 1
  
  kernel <- shapeKernel(c(3,3), type="disc")
  binROI = closing(binROI, kernel)
  binROI = opening(binROI, kernel)
  binROI = fillHull(binROI)
  
  ccROI = bwlabel(binROI)
  
  if (startslice){
    limiteTamanho = 1000
  }else{
    limiteTamanho = 200
  }
  rem = as.numeric(names(which(table(ccROI) < limiteTamanho)))
  for(i in rem){
    ccROI[ccROI == i] = 0
  }
  ccROI[ ccROI>1]= 1
  return (ccROI)
}

#### Calcula centroide das vertebras #####
calcula_centroides = function(mat){
  comp = bwlabel(mat)
  cent = c()
  for (i in names(table(comp)) ){
    v = as.numeric(i)
    m = comp
    m[m!=v] = 0
    posicoes = desindexa(m,which(m==v))
    centrox = mean(posicoes[1,])
    centroy = mean(posicoes[2,])
    cent = rbind(cent, c(centrox,centroy) )
  }
  cent = cent[-1,]
  return (cent)
}

#### Encontra o baco para a primeira fatia ####
acha_baco1 = function(m, vert){
  vertcent = calcula_centroides(vert)
  cent = calcula_centroides(m)
  if (length(cent) == 2){
    return(m)
  }
  ### Encontrar componente q minimiza soma das distancias as vertebras
  min_d = nrow(vertcent)*(2*nrow(m)^2)
  min_c = -1
  for (i in 1:nrow(cent)){
    d = 0
    x = cent[i,1]
    y = cent[i,2]
    for (j in 1:nrow(vertcent)){
      xv = vertcent[j,1]
      yv = vertcent[j,2]
      d = d+ (x-xv)^2 + (y-yv)^2
    }
    if (d < min_d){
      min_d = d
      min_c = i
    }
  }
  baco = bwlabel(m)
  baco[baco!=min_c] = 0
  baco[baco==min_c] = 1
  return (baco)
}

acha_baco = function(mat,baco1){
  cc = bwlabel(mat)
  max =0
  max_c = -1
  for (i in names(table(cc))[-1]){
    v = as.numeric(i)
    m = cc
    m[m!=v] = 0
    m[m==v] = 1
    comp = m * baco1
    if (length(table(comp))>1){
      pares = table(comp)[2]
    }else pares = 0
    if (pares > max){
      max = pares
      max_c = i
    }
  }
  m = cc
  m[m!=max_c] = 0
  m[m==max_c] = 1
  return(m)
}

#setwd("/Users/ludykong/MaChiron/MaChironGit")
home = getwd()# "/Users/ludykong/MaChiron/MaChironGit"

print("Comecando spleen segmentation!")

dir = "/Users/ludykong/GDrive/MaChiron/Exames/618999/DICOM/ARTERIAL/"
files = list.files(dir)
if (files[length(files)] == "Icon\r"){
  files = files[-length(files)]
}

# Encontrar o pixelSpacing
fname <-  paste(dir, files[1], sep="/")
dicom <- readDICOMFile(fname)
hdr = dicom$hdr

sl = which(hdr[,3]== "SliceLocation")
ipp = which(hdr[,3]== "ImagePositionPatient")
files[1]
hdr[ipp:sl,c(3,6)]

a = which(hdr[,3]== "PixelSpacing")
pixelSpacing = hdr[a,6]
b = unlist(strsplit(pixelSpacing, split = " "))
pixel_x = b[1]
pixel_y = b[2]
print(paste("pixel_x",pixel_x,"pixel_y",pixel_y))

#Recuperar ordem dos arquivos e fazer janela do baco
files = as.numeric(files)
files = files[order(files,decreasing = T)]
which(is.na(files))

if (length(unique(diff(files))) >1 ){
  print("Salto na altura das imagens!")
}

maximo = files[1]
relativos = maximo - files
COMECO_BACO = 700
FIM_BACO = 1500
BACO_BOM = 1000

bacos = relativos[relativos>COMECO_BACO & relativos<FIM_BACO]
baco1 = bacos[which(abs(bacos - BACO_BOM)==min(abs(bacos - BACO_BOM)) )]

baco1 = formatC(maximo-baco1, width = 4, format = "d", flag = "0")
bacos = formatC(maximo-bacos, width = 4, format = "d", flag = "0")

print(paste("N Bacos",length(bacos), "Baco bom:", which(bacos==baco1)))

bacos = c(baco1, bacos[-which(bacos==baco1)] )

startslice = TRUE
area_bacos = c()

file = bacos[1]
i=0
for (file in bacos){
    i=i+1
    print(paste(i,file))
    ori = le_dicom(dir,file)
#   EBImage::display(ori/max(ori))
    m = pre_normaliza(ori)
#   EBImage::display(m/max(m))
    vert = acha_brancos(m)
#   EBImage::display(vert/max(vert))
    m[vert==1]=0  
    if (startslice){
      jan = encontraROI(m,vert)
    }
    m = m[jan[1]:jan[2],jan[3]:jan[4]]
    vert = vert[jan[1]:jan[2],jan[3]:jan[4]]
    bin = filtra_componentes(m, startslice)
#   EBImage::display(bin)
    if (startslice){
      baco = acha_baco1(bin,vert)
      baco1 = baco
      startslice = FALSE
    } else{
      baco = acha_baco(bin,baco1)
    }
    EBImage::display(baco)
    tamanho_baco = table(baco)[2]
    #if (is.na(tamanho_baco)) break
    area_baco = tamanho_baco*as.numeric(pixel_x)*as.numeric(pixel_y)
    area_bacos = c(area_bacos,area_baco)
}
length(bacos)
area_bacos = unname(area_bacos)
print(area_bacos)

print("fim :)")

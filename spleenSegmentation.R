source("library.R")

##### Normaliza e pre processamento ######
pre_normaliza = function(m){
  bin = m > otsu(m,range = c(0,max(m)))
  binfill = fillHull(bin)
  binmaior = maior_componente(binfill)$matrix
  fabd = m*binmaior
  
  h = density(fabd[fabd>100])
  significativo = (h$y[-c(1,2)]>1e-4)
  
  i1 = h$x[which(diff(significativo)==1)]
  i2 = h$x[which(diff(significativo)==-1)]
  lo = i1
  hi = i2
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
  bin = m > otsu(m,range = c(0,max(m)))
  fB = m*bin
  
  muw = mean(fB[which(fB != 0)])
  sigw = sd(fB[which(fB != 0)])
  tw = muw + sigw 
  
  bin = fB>tw
  bin = fillHull(bin)
  kernel <- shapeKernel(c(3,3), type="disc")
  bin_pos = opening(bin, kernel)
  bin_pos = closing(bin_pos, kernel)
  
  # Remover componentes pequenos da mascara
  limiteTamanho = 20
  cc = bwlabel(bin_pos)  
  
  rem = as.numeric(names(which(table(cc) < limiteTamanho)))
  for(i in rem){
    cc[cc == i] = 0
  }
  cc[ cc >1]= 1
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
  i=3
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

dir = "/Users/ludykong/GDrive/MaChiron/Exames/Cisto/DICOM/VOLUME ABD SC 1.0  B20f"
print(paste("Pasta:",dir))

files = list.files(dir)
if (files[length(files)] == "Icon\r")
  files = files[-length(files)]

fname <-  paste(dir, files[1], sep="/")
dicom <- readDICOMFile(fname)
hdr = dicom$hdr
a = which(hdr[,3]== "PixelSpacing")
pixelSpacing = hdr[a,6]
b = unlist(strsplit(pixelSpacing, split = " "))
pixel_x = b[1]
pixel_y = b[2]
print(paste("pixel_x",pixel_x,"pixel_y",pixel_y))

startslice = TRUE

area_bacos = c()

file = "0006-65662095"
file = files[1]
for (file in files[10:14]){
    ori = le_dicom(dir,file)
    display(ori/max(ori))
    m = pre_normaliza(ori)
    vert = acha_brancos(m)
    m[vert==1]=0  
    if (startslice){
      jan = encontraROI(m,vert)
    }
    m = m[jan[1]:jan[2],jan[3]:jan[4]]
    vert = vert[jan[1]:jan[2],jan[3]:jan[4]]
    m = filtra_componentes(m, startslice)

    if (startslice){
      baco = acha_baco1(m,vert)
      baco1 = baco
      startslice = FALSE
    } else{
      baco = acha_baco(m,baco1)
    }
    display(baco)
    tamanho_baco = table(baco)[2]
    area_baco = tamanho_baco*as.numeric(pixel_x)*as.numeric(pixel_y)
    area_bacos = c(area_bacos,area_baco)
}

print(area_bacos)

print("fim :)")

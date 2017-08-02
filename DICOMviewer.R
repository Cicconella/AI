source("library.R")
print("Comecando spleen segmentation!")

dir = "/Users/ludykong/GDrive/MaChiron/Exames/1775933/DICOM/ARTERIAL/"
files = list.files(dir)
if (files[length(files)] == "Icon\r"){
  files = files[-length(files)]
}
fname <-  paste(dir, files[1], sep="/")
dicom <- readDICOMFile(fname)
hdr = dicom$hdr
a = which(hdr[,3]== "PixelSpacing")
pixelSpacing = hdr[a,6]
b = unlist(strsplit(pixelSpacing, split = " "))
pixel_x = b[1]
pixel_y = b[2]
print(paste("pixel_x",pixel_x,"pixel_y",pixel_y))
a = which(hdr[,3]== "SliceThickness")
st = as.numeric(hdr[a,6])/2
print(paste("Slice T",st ) )
k = as.numeric(st)/as.numeric(pixel_x)
x =512
z = 267
m = dicom$img
dim(m)
length(files)
imagem3D = array(0,c(nrow(m),ncol(m),length(files)))

i=0
for (file in files){
  i=i+1
  print(paste(i,file) )
  fname <-  paste(dir, file, sep="/")
  dicom <- readDICOMFile(fname)
  m = dicom$img
  imagem3D[,,i] = m
}

D = dim(imagem3D)
#Axial
axial = array(0,c(D[2],D[1],D[3]))
for (corte in 1:D[3]){
  for (li in 1:D[1]){
    axial[,1+D[1]-li,corte] = imagem3D[li,,corte]
  }
}
EBImage::display(axial/max(axial),tile = "Axial")
#EBImage::display(imagem3D/max(imagem3D),method = "raster", title="All images",all=T,nx=14)

#Frontal
frontal = array(0,c(D[2],D[3],D[1]))
for (corte in 1:D[2]){
  for (li in 1:D[3]){
    frontal[,1+D[3]-li,corte] = imagem3D[corte,D[1]:1,li]
  }
}
EBImage::display(frontal/max(frontal))

#Sagital
sagital = array(0,c(D[1],D[3],D[2]))
for (corte in 1:D[2]){
  for (col in 1:D[3]){
    sagital[,1+D[3]-col,corte] = imagem3D[,corte,col]
  }
}
EBImage::display(sagital/max(sagital))


png("teste.png",width = 512, height=338)
image(t(ma), col=grey(0:64/64), axes=FALSE, xlab="", ylab="")
dev.off()
#png("teste.png",width = 512, height=338)
#image(t(ma), col=grey(0:64/64), axes=FALSE, xlab="", ylab="")
#dev.off()


plot(ma)
m = array(c(0,0.3,0.6,1),c(10,10,2,3))
m[1,1,1,1]
EBImage::display(x = array(ma/max(ma),ma/max(ma)),nx=2)
EBImage::display(m/max(m),all=T)

#install.packages("oro.nifti")
source("library.R")
library(oro.nifti)

#fname = "/Users/ludykong/Downloads/segmentation-0.nii"
fname = "/Users/ludykong/Desktop/volume-0.nii"
a = matrix(c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4),nrow=4,byrow = T)
write.table(a,"/Users/ludykong/Desktop/Tabela.txt",col.names=F,row.names=F)

NII <- readNIfTI(fname)
writeNIfTI(NII,paste0(path,"/NIFTISPHER"))
pixdim(NII)
oro.nifti::dim_(NII)
M <- oro.nifti::img_data(NII)

path = "/Users/ludykong/Desktop/TabelaNIFTI"
if(dir.exists(path)){
  print(paste("O diretorio", path,  "ja existe."))
  unlink(path,recursive = T)  
}else{
  dir.create(path)
}
ori = M[,,50]
plota_imagem(ori)

for (i in 1:dim(M)[3]){
  write.table(M[,,i],paste0(path,'/slice-',i),,col.names=F,row.names=F)
}

dir = "/Users/ludykong/Desktop/ARTERIAL/"
files = list.files(dir)
files[length(files)]
unlink(paste0(dir,"Icon\r"))
if (files[length(files)] == "Icon\r"){
  files = files[-length(files)]
}
?readDICOM
D = readDICOM(dir, exclude = "Icon\r", verbose=TRUE)
file = readDICOMFile(paste0(dir,files[1]))
M = D$img

dcmSphere <- readDICOM(system.file("sphere3", package="oro.dicom"), verbose=TRUE)
writeHeader(dcmSphere$hdr,"headder")
M = dcmSphere$img

ori = M$`/Users/ludykong/GDrive/MaChiron/Exames/1775933/DICOM/ARTERIAL//11245.dcm`
dim(ori)
plota_imagem(ori)
ni = dicom2nifti(D)
M = img_data(ni)
ni2 = oro.nifti::`pixdim<-`(ni,c(-1,1, 1,1, 1,1, 1,1))
writeNIfTI(ni,"/Users/ludykong/Desktop/NIFIII")
?oro.dicom::dicom2nifti
convert.datatype()

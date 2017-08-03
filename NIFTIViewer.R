#install.packages("oro.nifti")
library(oro.nifti)
library(EBImage)
#fname = "/Users/ludykong/Downloads/segmentation-0.nii"
fname = "/Users/ludykong/Desktop/volume-0.nii"
a = matrix(c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4),nrow=4,byrow = T)
write.table(a,"/Users/ludykong/Desktop/Tabela.txt",col.names=F,row.names=F)

NII <- readNIfTI(fname)
pixdim(NII)
oro.nifti::dim_(NII)
M <- oro.nifti::img_data(NII)
m = M[1:4,1:4,1:6]

path = "/Users/ludykong/Desktop/TabelaNIFTI"
if(dir.exists(path)){
  print(paste("O diretorio", path,  "ja existe."))
}else{
  dir.create(path)
}
ori = M[,,50]
EBImage::display(ori/max(ori))
for (i in 1:dim(M)[3]){
  write.table(M[,,i],)
}
#dir.exists(path)
#unlink(pathImagens,recursive=T)

image(NII,plot.type="single",z=60)


#### Libraries and functions ####
library(oro.dicom)
library(imager)
library(mmand)
library("e1071")
#library(cluster) Old FCM Clustering

#### Abrindo o arquivo ####

print("Iniciando Script Liver Detect")
base_dir = "/Users/ludykong/MaChiron/"
file_dir = "Data/52490002/"
arquivos <- list.files(paste(base_dir,file_dir, sep=""))
print(arquivos)

#print(paste("Diretorio:",base_dir))
#print(paste("Fatia alvo:",arquivo))

# fotos interessantes
#65647276 + 3 proximos 
#65646836
#65646946
#65652505

# Outros arquivos #
#65643294 
#65643316
#65643338
#65644086
#65649322

#64218107
#arquivos[9:20]
for (arquivo in arquivos[121:140]){
  fname <- paste(base_dir,file_dir,arquivo,sep="")
  print(fname)
  dicom <- readDICOMFile(fname)
  #names(dicom)
  #head(dicom$hdr)
  dim(dicom$img)
  png(paste(fname,".png",sep=""),width=1024,height=1024)
  image(t(dicom$img), col=grey(0:64/64), axes=FALSE, xlab="", ylab="")
  dev.off()
}
print( paste(c("Dicom eh imagem com dimensoes:",dim(dicom$img)),collapse=" ") )

#### Histograma Original e Normalizacao ####
#Valores Maximos e Minimos
#min(apply(a, 1, min))
#max(apply(a, 1, max))
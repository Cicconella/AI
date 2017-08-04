
library(oro.dicom) # Leitura do arquivo DICOM

dir = "/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/"

filename = "MaChiron/Exames/1775933/DICOM/ARTERIAL/12865.dcm"
fname <-  paste(dir, filename, sep="")

abdo <- readDICOMFile(fname)

plota_imagem(abdo)
hdr = abdo$hdr
ori = abdo$img
#plot_image(ori, "Original")
plota_imagem(ori)

hist(ori)

head(hdr[,c(3,6)])

hdr[c(which(hdr[,3]=="RescaleIntercept"),which(hdr[,3]=="RescaleSlope")),6]

ori = ori* as.numeric(hdr[which(hdr[,3]=="RescaleSlope"),6])+
  as.numeric(hdr[which(hdr[,3]=="RescaleIntercept"),6])

hist(ori)

ori[which(ori<(-100))]=-1000
ori[which(ori>400)]=-1000

plota_imagem(ori)

ori2= equalize(ori, range = c(min(ori), max(ori)), levels = 256)

plota_imagem(ori2)

hist(ori2[ori2!=(-1000)])

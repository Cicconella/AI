
install.packages("oro.dicom")
library(oro.dicom)

### Script para adicionar um novo exame ao banco de dados

### Dado que getwd() = "/home/cicconella/AI", meu banco vai ficar aqui:

database = paste(getwd(), "/../GugolDraive/MaChiron/Exames", sep="")

list.files(database)

exame = "Downloads"

path = paste0(database,"/", exame)

if(dir.exists(path)){
  print(paste("O diretorio", path,  "ja existe."))
}else{
  dir.create(path)
  pathDICOM = paste0(path, "/DICOM")
  pathImagens = paste0(path, "/Imagens")
  
  dir.create(pathDICOM)
  dir.create(pathImagens)
}

# Diretorio de uma pasta de exame contendo muitos DICOMs
dir = paste(getwd(), "/../", exame,"/DICOM/", sep="")

#Lista os DICOM
nomes = list.files(dir)

#nomes = read.table(paste(dir,"nome",sep=""))
#dim(nomes)
#nomes = nomes[-dim(nomes)[1],1]
#nomes = as.character(nomes[,1])
#head(nomes)
#nomes[807:818]
#nomes =  nomes[-c(807:818)]

paciente = NA
estudo = NA
tipo = rep(NA, length(nomes))
slice = rep(NA, length(nomes))
fase = rep(NA, length(nomes))

i=1
for(i in 1:length(nomes)){
  fname <- paste(dir, nomes[i], sep="")
  
  abdo <- readDICOMFile(fname)
  header = abdo$hdr[,c(3,6)]

  if(i==1){
    paciente = header[which(header[,1]== "PatientID"),2]
    estudo = header[which(header[,1]== "StudyDescription"),2]
    espessura = header[which(header[,1]== "SliceThickness"),2]
  }
  
  tipo[i] = header[which(header[,1]== "ImageType"),2]
  
  a = which(abdo$hdr[,3]== "SliceLocation")
   if(length(a)!=0){
     slice[i] = (abdo$hdr[a,6])
   }else{
     slice[i] = NA
   }
  
  a = which(abdo$hdr[,3]== "SeriesDescription")
  if(length(a)!=0){
    fase[i] = (abdo$hdr[a,6])
  }else{
    fase[i] = NA
  }
  
   png(paste(pathImagens,"/",nomes[i], ".png", sep=""),width = 500,height = 500)
   image(t(abdo$img), col=grey(0:64/64), axes=FALSE, xlab="", ylab="")
   dev.off()
  
  print(i)
  
}


info = cbind(nomes, tipo, fase, slice)
head(info)

write.table(rbind(c("PAX: ", paciente), c("Tipo de Exame: ", estudo),c("Espessura: ", espessura)), paste(path,"/infosExame", sep=""), 
            row.names = F, quote = F, col.names = F)
write.table(cbind(nomes, tipo, fase, slice), paste(path,"/infosImagem", sep=""), row.names = F)

sc = grep("SC", fase)
art = grep("ART", fase)
eq = grep("EX", fase)
pont = seq(1,length(fase))[-c(sc, art,eq)]


nfases  = unique(fase)
for (f in nfases){
  dir.create(paste0(pathDICOM,"/",f))
} 


for (i in 1:length(nomes)){
  dicom = nomes[i]
  f = fase[i]
  s = formatC(as.numeric(slice[i])*10, width = 4, format = "d", flag = "0")
  file.rename(paste0(dir,"/",dicom),paste0(pathDICOM,"/",f,"/",s,"-",dicom) )
}

unlink(dir,recursive=T)
dir.exists(dir)

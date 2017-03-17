library(oro.dicom)

dir = "/media/cicconella/8AA6013CA6012A71/Documents and Settings/Nina/Google Drive/MaChiron/Exames/HNF grande/"

nomes = read.table(paste(dir,"nome",sep=""))
dim(nomes)
nomes = nomes[-dim(nomes)[1],1]
nomes= as.character(nomes)
head(nomes)

paciente = rep(NA, length(nomes))
tipo = rep(NA, length(nomes))
estudo = rep(NA, length(nomes))
slice = rep(NA, length(nomes))
fase = rep(NA, length(nomes))

for(i in 1:length(nomes)){
  fname <- paste(dir, "DICOM/", nomes[i], sep="")
  
  # if(i==809){
  #   paciente[i] = NA
  #   tipo[i] = NA
  #   estudo[i] = NA
  #   slice[i] = NA
  #   fase[i] = NA
  #   next
  # }
  
  abdo <- readDICOMFile(fname)
  
  abdo$hdr[,c(3,6)]

  a = which(abdo$hdr[,3]== "PatientID")
  paciente[i] = (abdo$hdr[a,6])

  # abdo$hdr[abdo$hdr$name == "SliceThickness",] #64
  
  a = which(abdo$hdr[,3]== "ImageType")
  tipo[i] = (abdo$hdr[a,6])
  
  a = which(abdo$hdr[,3]== "StudyDescription")
  estudo[i] = (abdo$hdr[a,6])
  
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
  
#   png(paste(dir,"Imagens/",nomes[i], ".png", sep=""))
#   image(t(abdo$img), col=grey(0:64/64), axes=FALSE, xlab="", ylab="")
#   dev.off()
  
  print(i)
  
}

write.table(cbind(paciente, estudo, tipo, slice, fase, nomes), paste(dir,"infos", sep=""), row.names = F)





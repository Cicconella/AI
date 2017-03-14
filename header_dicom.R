library(oro.dicom)

dir = "/media/cicconella/8AA6013CA6012A71/Documents and Settings/Nina/Dropbox/Machiron Dicoms/Hemangioma grande lobo esquerdo/iSiteExport/DICOM/"

nomes = read.table(paste(dir,"nome",sep=""))
dim(nomes)
nomes = nomes[-dim(nomes)[1],1]
head(nomes)

pacienteNome = rep(NA, length(nomes))
paciente = rep(NA, length(nomes))
slice = rep(NA, length(nomes))
fase = rep(NA, length(nomes))

for(i in 1:length(nomes)){
  fname <- paste(dir, nomes[i], sep="")
  
  abdo <- readDICOMFile(fname)
  abdo$hdr[,c(3,6)]

  
  a = which(abdo$hdr[,3]== "PatientsName")
  paciente[i] = (abdo$hdr[a,6])
  
  a = which(abdo$hdr[,3]== "PatientID")
  paciente[i] = (abdo$hdr[a,6])

   a = which(abdo$hdr[,3]== "SliceLocation")
   if(length(a)!=0){
     slice[i] = (abdo$hdr[a,6])
   }else{
     slice[i] = NA
   }
  
  
   a = which(abdo$hdr[,3]== "SeriesDescription")
   fase[i] = (abdo$hdr[a,6])
  
  
   png(paste(dir,"Imagens/",nomes[i], ".png", sep=""))
   image(t(abdo$img), col=grey(0:64/64), axes=FALSE, xlab="", ylab="")
   dev.off()
  
  print(i)
  
}

paciente

paciente[1:100]
slice[1:100]
fase[1:100]

write.table(cbind(nomes, paciente, slice, fase), paste(dir,"infos", sep=""), row.names = F)

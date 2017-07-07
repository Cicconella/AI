dir = "/Users/ludykong/GDrive/MaChiron/Exames/124594/DICOM/ARTERIAL/"
files = list.files(dir)
for (file in files){
  x =unlist(strsplit(file,"-"))
  if (length(x)==3){
    a = "-"
    novo_nome = paste(a,x[2],sep='')
  }else{
    novo_nome = x[1]
  }
  file.rename(paste0(dir,file),paste0(dir,novo_nome) )
}

file = files[1]
for (file in files){
  x =unlist(strsplit(file,".dcm"))
  novo_nome = x[1]
  file.rename(paste0(dir,file),paste0(dir,novo_nome) )
}

file
files = list.files(dir)
if (files[length(files)] == "Icon\r"){
  files = files[-length(files)]
}
file = files[1]
i=0
for (file in files){
  i=i+1
  print(paste(i,length(files), file))
  fname <-  paste0(dir, file)
  dicom <- readDICOMFile(fname)
  hdr = dicom$hdr

  ipp = which(hdr[,3]== "ImagePositionPatient")
  x = hdr[ipp,c(6)]
  b = unlist(strsplit(x, split = " "))
  z = as.numeric(b[3])
  z4 = formatC(round(10*z), width = 4, format = "d", flag = "0")
  file.rename(paste0(dir,file),paste0(dir,z4) )
}

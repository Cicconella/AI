##### Bibliotecas #####
source("library.R")

##### Nossos dados ###### 
#dir da ana
dir = "/media/cicconella/01D2FE834EA51BE0/Documents and Settings/Nina/Google Drive/"

#dir do luis
#dir = "/Users/ludykong/GDrive/"

filename = "MaChiron/Exames/Teste LITS/volume-0.nii"

fname <-  paste(dir, filename, sep="")

abdo <- readNIfTI(fname)
abdo
image(abdo,plot.type="single",z=60)
#hist(ori[ori>10], nc=1000,main = "Histograma Original")


m = img_data(abdo)
m = m[,,60]

hist(m)
teste = density(m)
teste = as.vector(m)
teste = teste[teste>-2000]

hist(teste)

m = m + abs(min(m))

hist(m)


plota_imagem(m)

norm = pre_normaliza(m)
plota_imagem(norm)

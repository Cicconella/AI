#install.packages("oro.nifti")
library(oro.nifti)
library(EBImage)
fname = "/Users/ludykong/Downloads/segmentation-0.nii"
fname = "/Users/ludykong/Downloads/volume-0.nii"
abdo2 <- readNIfTI(fname)
image(abdo2,plot.type="single",z=60)
image(abdo,plot.type="single",z=60)
#plot_image(ori, "Original")
EBImage::display(ori/max(ori))

setwd("Desktop/Dicom")
library("oro.dicom")
library("dcemriS4")

## read dicom
setwd("/home/pang/data/SCP_involve/Dicom/XA_Untersuchungen_test")
library(oro.dicom)
dcmImages <- readDICOM("Abels Josef 13.10.1955 - XA 24.05.2019 D2000_19", verbose = TRUE,recursive = FALSE, exclude = "sql")
dcm.info <- dicomTable(dcmImages$hdr)
length(names(dcm.info))
unique(dcm.info["0008-0008-ImageType"])
unique(dcm.info["0018-1110-DistanceSourceToDetector"])
unique(dcm.info["0018-1130-TableHeight"])
image(t(dcmImages$img[[1]]), col = grey(0:64/64), axes = FALSE,
      xlab = "", ylab = "")
extractHeader(dcm$hdr, "Manufacturer", numeric=FALSE)
image(t(dcm$img), col=grey(0:64/64), axes=FALSE, xlab="", ylab="",
      main="Example image from DICOM file")
system.time(for (stack in unique(substring(rownames(dcm.info), 15, 20))) {
  if (substring(stack, 1, 2) == "CT") {
    print(stack)
    index <- which(substring(rownames(dcm.info), 15, 20) == stack)
    dcm.stack <- list(hdr = dcmImages$hdr[index], img = dcmImages$img[index])
    dcm.nifti <- dicom2nifti(dcm.stack, DIM = 3, descrip = c("Manufacturer",
                                                             "ManufacturersModelName"))
    print(dcm.nifti)
    writeNIfTI(dcm.nifti, paste0("NIfTI/", stack))
  }
})
orthographic(dcm.nifti, col.crosshairs = "green")
image(dcm.nifti)



dcm<-readDICOM("test")
dcm<-readDICOM("PAT47838948",debug=TRUE)
dcmImages2 <- readDICOM("test", verbose = TRUE,recursive = FALSE, exclude = "sql")
dcmImages <- readDICOM("PAT47838948/test", verbose = TRUE,recursive = FALSE, exclude = "sql")
dcm.info <- dicomTable(dcmImages$hdr)
length(names(dcm.info))
unique(dcm.info["0008-0008-ImageType"])
unique(dcm.info["0018-1110-DistanceSourceToDetector"])
unique(dcm.info["0018-1130-TableHeight"])
image(t(dcmImages$img[[1]]), col = grey(0:64/64), axes = FALSE,
      xlab = "", ylab = "")
extractHeader(dcm$hdr, "Manufacturer", numeric=FALSE)
#image(t(dcm$img), col=grey(0:64/64), axes=FALSE, xlab="", ylab="",
      main="Example image from DICOM file")
for (stack in unique(substring(rownames(dcm.info), 15, 20))) {
  if (substring(stack, 1, 2) == "CT") {
    print(stack)
    index <- which(substring(rownames(dcm.info), 15, 20) == stack)
    dcm.stack <- list(hdr = dcmImages$hdr[index], img = dcmImages$img[index])
    dcm.nifti <- dicom2nifti(dcm.stack, DIM = 3, descrip = c("Manufacturer",
                                                             "ManufacturersModelName"))
    print(dcm.nifti)
    writeNIfTI(dcm.nifti, paste0("NIfTI/", stack))
  }
}
orthographic(dcm.nifti, col.crosshairs = "green")
image(dcm.nifti)

for (i in 1:33){
  png(paste("PAT47838948/cordcm/image",i,".png",sep=""))
  image(t(dcmImages$img[[1]][,,i]), col = grey(0:64/64), axes = FALSE, xlab = "", ylab = "")
  dev.off()
}

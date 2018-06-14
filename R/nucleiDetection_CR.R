library(CRImage)
library(png)
files <- list.files("/home/aykan/layer0patches/origImagesHE", full.names = TRUE, recursive = TRUE)
count = 0
for (f in files) {
  print(f)
  count = count + 1
  segVals = segmentImage(filename=f, maxShape = 800, minShape = 40, failureRegion = 2000, numWindows = 4)

  write.table(segVals[[3]], file = paste("NucleiR_HE", toString(count), ".txt"))
  writePNG(segVals[[2]], target = paste("Mask_HE", toString(count), ".png"))
  print(count)
}

files <- list.files("/home/aykan/layer0patches/origImagesHE", full.names = TRUE, recursive = TRUE)
count = 0
for (f in files) {
  print(f)
  count = count + 1
  segVals = segmentImage(filename=f, maxShape = 800, minShape = 100, failureRegion = 2000, numWindows = 4)

  write.table(segVals[[3]], file = paste("NucleiR100_HE", toString(count), ".txt"))
  writePNG(segVals[[2]], target = paste("Mask100_HE", toString(count), ".png"))
  print(count)
}

files <- list.files("/home/aykan/layer0patches/origImagesTP53", full.names = TRUE, recursive = TRUE)
count = 0
for (f in files) {
  print(f)
  count = count + 1
  segVals = segmentImage(filename=f, maxShape = 800, minShape = 100, failureRegion = 2000, numWindows = 4)

  write.table(segVals[[3]], file = paste("NucleiR100_TP", toString(count), ".txt"))
  writePNG(segVals[[2]], target = paste("Mask100_TP", toString(count), ".png"))
  print(count)
}

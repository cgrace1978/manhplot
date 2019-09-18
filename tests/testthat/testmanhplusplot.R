context("Running the manhplusplot on (small) dummy data")
library(manhplot)

test_that("Run the manhplusplot function with default params", {
  infile<-test_path("cad.add.160614_manhformat_chr1.txt.gz")
  configfile<-test_path("config.txt")
  snpfile<-test_path("5cad.add.160614.variants_chr1.txt")

  manhplusplot(infile = infile, outfile = file.path(tempdir(), "testpdf"),configfile = configfile, snpfile = snpfile)
})

test_that("Run the manhplusplot function with output as tiff file", {
  infile<-test_path("cad.add.160614_manhformat_chr1.txt.gz")
  configfile<-test_path("config.txt")
  snpfile<-test_path("5cad.add.160614.variants_chr1.txt")
  
  manhplusplot(infile = infile, outfile = file.path(tempdir(), "testtiff"),configfile = configfile, snpfile = snpfile, drawastiff = T)
})

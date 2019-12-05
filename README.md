# MANHATTAN++

MANHATTAN++ is software to generate a transposed manhattan heatmap, implemented in R.

## Getting Started

You need to install the latest version of [R](https://www.r-project.org/) The R package can be run on Windows and Linux, you must specify paths to the filenames you are using as input and output.

To install the software from the GIT repository:
```
install.packages("devtools")
library(devtools)

install_github("cgrace1978/manhplot", dependencies = T, force = T)
```

The following command will run the plot with default data in the package. The pdf (test.pdf) will be created in the current working directory in R (This can be viewed using the getwd() command):
```
library(manhplot)

infile<-system.file("extdata","cad.add.160614_manhformat.txt.gz",package = "manhplot")
configfile<-system.file("extdata","config.txt", package = "manhplot")
snpfile<-system.file("extdata","56cad.add.160614.variants.txt", package = "manhplot")

## Run manhattan++ with the default paramaters and files included in the package
manhplusplot(infile = infile,outfile = "test", configfile = configfile, snpfile = snpfile)
```
For more information on using manhattan++ please visit the [Manhattan++ wiki](https://github.com/cgrace1978/manhplot/wiki/home)



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

The following command will run the plot with default data in the package:
```
library(manhplot)

infile<-system.file("extdata","cad.add.160614_manhformat.txt.gz",package = "manhplot")
configfile<-system.file("extdata","config.txt", package = "manhplot")
snpfile<-system.file("extdata","56cad.add.160614.variants.txt", package = "manhplot")

## Run manhattan++ with the default paramaters and files included in the package
manhplusplot(infile = infile,outfile = "test", configfile = configfile, snpfile = snpfile)
```
## Input Files

In order to generate the plot, three files are required (with headers in the format described)

### Genome Wide Association Study (GWAS) or Meta-analysis summary statistics

The summary statistics, this file should have the following columns:

1. **chr** - chromosome (must be numeric)
2. **pos** - position
3. **pvalue** - p-value (please **do not** log-transform)
4. **maf** - Minor Allele Frequency
5. **conseq** - Flag of whether variant has HIGH functional consequence (using VEP)

The variable *infile* should be modfied to point to this file (with the file structure of the operating system used). The files should be numerically sorted (chr:pos).
```
chr pos     pvalue      maf       conseq
1   751756  0.4528019   0.158264  0
1   752566  0.7394597   0.763018  0
1   752721  0.8462652   0.740969  0
1   752894  0.7750657   0.744287  0
```
**HINT:** The software's default is to displays MAF data when a MAF <5%. If you don't want to use this display function set the MAF function argument to 0.

**HINT:** The software displays annotation data when a consequence is 1. If you don't want to use this display function, give all SNPs a consequence of 0

**HINT:** You can use different column names for the GWAS files by using the command line arguments:
```
## Use different GWAS file column names using: chrname, posname, pvalname, frqname and conseqname
manhplot(infile = infile,outfile = "test", configfile = configfile, snpfile = snpfile,
         chrname = "chr", posname = "pos", pvalname = "pvalue", frqname = "maf",
         conseqname = "conseq",)
```

### Loci information (variants)

Variants of interest to annotate the plot with. 

1. **markername** - The name of the variant
2. **NearestGene** - The gene name associated with the variant
3. **chr** - chromosome
4. **pos** - position
5. **eaf** - effect allele frequency
6. **OR** - Odds Ratio
7. **Pvalue** - P-value
8. **novel** - Flag for whether the variant is novel (a new loci)

The variable *snpfile* should be modfied to point to this file (with the file structure of the operating system used). The novel flag indicates whether the loci is a new finding in the GWAS being reported. When the *showGenes* flag is set to TRUE, Novel loci will be displayed in the table and the others on the heatmap, as bubble labels. When *showGenes* is false, all loci are displayed in the table.

```
markername  NearestGene chr     pos         eaf       OR      Pvalue    novel
rs11206510  PCSK9       1       55496039    0.847627  1.08    2.34E-08  FALSE
rs9970807   PPAP2B      1       56965664    0.915097  1.13    5E-14     FALSE
rs7528419   SORT1       1       109817192   0.78582   1.12    1.97E-23  FALSE
```

### Configuration file

To define the legend and type of heatmap cells.

1. **min.count** - Lower bounds of number of variants in this cell
2. **max.count** - Upper bounds of number of variants in this cell.
3. **maf** - Should variants with minor allele frequencies below the threshold be detected by this cell?
4. **conseq** - Should any HIGH consequence variant be detected?
5. **col** - The colour of the heatmap cells
6. **idx** - Index for labeling on the heatmap
7. **type** - The type of cell: val or oddchr / evenchr
8. **report** - Should these cells be reported on the heatmap with bubble label?

The variable *configfile* should be modfied to point to this file (with the file structure of the operating system used)
```
#### CONFIG file for use with MANH++ - Do not modify the first 10 lines							
## min.count: The lower cell count threshold to accept this config							
## max.count: The upper cell count threshold to accept this config							
## maf: Is MAF detection active for this config - is there any variants within a cell with MAF < threshold?	
## conseq: Is HIGH impact consequence active? Are there any variants with HIGH impact consequence in the cell?		
## col: The colour which cells for this config							
## idx: the index to use for this cell in the heatmap - MUST BE CONSECUTIVE FROM START TO END - STARTING AT 1
## type: val - an config entry, oddchr - the odd chromosome, evenchr - the even chromosome"		
## report: Are these annotations labeled on the heatmap							
#####		
min.count   max.count   maf     conseq    col         idx     type    report
1           2           FALSE   FALSE     black       1       val     FALSE
1           2           FALSE   TRUE      light pink  2       val     TRUE
1           2           TRUE    FALSE     green       3       val     FALSE
```

Also see the example file (config.txt) for more information.

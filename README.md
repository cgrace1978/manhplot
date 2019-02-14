# MANHATTAN++ v1.0

MANHATTAN++ is software to generate a transposed manhattan heatmap, implemented in R.

## Getting Started

you need to install [R v3.4.3](https://www.r-project.org/) and the packages **ggplot2**, **reshape2**, **ggrepel** and **gridExtra** and associated dependencies. The R script can be run on Windows and Linux, you must specify paths to the filenames you are using (or use the default filenames)

To install the software you need to download this GIT repository, this can be done on the command line using the command:
```
git clone https://github.com/cgrace1978/manhplot
```
To update the local version of the repository use the following command:
```
cd manhplot
git pull
```
Alternatively it is possible to download the repository as a zip using the clone / download on the GIT UI.

Once the repository is downloaded, install the required packages in R and set the file paths to your input files, and the script should just work with the command:
```
Rscript manhattan.heatmap.v1.R
```
**Note:** If you are using the default GWAS file (cad.add.160614_manhformatv3.txt.gz) it must be unzipped prior to running the R script.
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
**HINT:** The software displays MAF data when a MAF <5%. If you don't want to use this display function, give all SNPs a MAF > 5%

**HINT:** The software displays annotation data when a consequence is 1. If you don't want to use this display function, give all SNPs a consequence of 0

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

## Configuration options (within script)

### Flags
The script has the following flags:

1. **showgenes** - Set to TRUE to show labels for known genes, rather than cells of interest
2. **showrsids** - Set to TRUE to show rsids rather than genes on labels
3. **drawastiff** - Set to TRUE to draw the plot as a TIFF file (default is PDF)

### Variables
The following variables can be customized:

1. **pos.split** - The length of base pair regions for each cell - default 3E6
2. **pval.split** - The length of log10(pvalue) for each cell - default 0.125
3. **max.pval** - the max pval to show on the plot - default 20
4. **pval.units** - The -log10(pval) breaks to display on the y axis
5. **textsize** - The size of text used on labels
6. **GWS** - Genome wise significant threshold to use in the plot - default 5E-8
7. **FDR** - The FDR 5% threshold to use on the plot - default 6.28E-5, Below this p-value cells are ignored and greyed out
8. **MAF** -  The MAF threshold for use with the config MAF flag.

### File variables
The following variables point to input / output files (modify to operating system file format)

1. **infile** - File containing the GWAS summary statistics
2. **outfile** - Prefix of the output file (PNG or TIFF) (default date stamp)
3. **snpfile** - File containing Loci information to be displayed on the figure.
4. **configfile** - Configuration file - determines colour scheme and associated variant bins for the cells of the heatmap.
5. **debugfile** -File with debug / logging, activated using the *debugflag* flag. Writes the contents of each cell in the heatmap

## R version 3.4.2 or higher required

library(reshape2)
library(ggplot2)
library(ggrepel)
require(gridExtra)

## files by default are located in respect to the working directory for this project
## use getwd and setwd to change this.
## input gwas file needs to be numerically sorted with columns:
## chr, pos, Pvalue, allele frequency, VEP annotation
## VEP annotation coded as 1 if present at variant, 0 is not present
## allele frequency is reported if <- MAF variable below.
infile<-"cad.add.160614_manhformatv3.txt" ## ukbb.full.txt cad.add.160614_manhformatv2.txt

## name of the PDF to generate - currently using a date stamp - modify if necessary.
outfile<-paste("manh-",format(Sys.time(), "%d-%m-%y.%H-%M-%S"),sep="")

## name of the config file
configfile<-"config.txt"

## snpfile needs contains the following columns:
## markername, eaf, OR, NearestGene, chr, pos, Pvalue
## snp list should be length less than 130 to be displayed optimally
snpfile<-"56cad.add.160614.variants.txt" ## snp.infov2.txt

## name of log file to generate if debugflag is set.
debugfile<-paste("manh_",format(Sys.time(), "%d-%m-%y.%H-%M-%S"),".log",sep="")

## FLAGS
showgenes<-F ## show labels for known genes, rather than cells of interest
showrsids<-F ## show rsids on labels
drawastiff<-F ## Draw the plot as a TIFF file (default is PDF)

## parameters for drawing the manhattan heatmap
pos.split<-3E6 ## length of base pair regions for each cell - default 3E6
pval.split<-0.125 ## length of lo10(pvalue) for each cell - default 0.125
max.pval<-20 ## the max pval to show on the plot - default 20
pval.units<-5 ## units to display on the y axis
textsize<-2 ## size of text used on labels

## significance thresholds to display
GWS<-5E-8
FDR<-1E-3 ## All cells below this threshold are ignored and greyed out

MAF<-0.05 ## The MAF threshold for use with the config MAF flag.

#### THERE IS NO NEED TO READ THE SCRIPT BEYOND THIS POINT

rebuild<-F ## set to false to retain the current matrix
debugflag<-F ## Turn logging on / off

### Assert statement
## condition to test
## message if the condition fails.
waitifnot <- function(cond, mess) {
  if (!cond) { ## check that condition is fulfilled
    message(mess)
    message(deparse(substitute(cond)), " is not TRUE")  
    while (TRUE) {}
  }
}

### function returns y index on the heatmap that maps to a specific -log10(p-value)
## log - log10 pval to convert to index
log10.index<-function(log){
  idx<-(log / pval.split) + 0.5
  return(idx)
}

### function returns y index on the heatmap that maps to a specified p-value 
## p.val - pval to convert to index
p.val.index<-function(p.val){
  gws<--1*log10(p.val) ## convert the p-value to -log10
  idx<-log10.index(gws) ## call the log10.index function.
  return(idx)
}

p.val.cell<-function(p.val,pval.chunks){
  gws<--1*log10(p.val) ## convert the p-value to -log10
  idx<-log10.cell(gws) ## call the log10.index function.
  return(idx)
}

log10.cell<-function(log, pval.chunks){
  idx<-pvals.cells.index$id[log >= pvals.cells.index$LP & log < pvals.cells.index$UP]
  return (idx)
}

### Find the exact cell which the position is within on the heatmap
## chr - chr for index
## position - position for index
## chr.chunks - datastructure containing cell locations for chromosomes and positions.
chrpos.cell<-function(chr, position, chr.chunks){
  slice<-chr.chunks[chr.chunks$chr==chr & chr.chunks$s < position & chr.chunks$e > position,]
  
  return (slice$posid)
}

if(rebuild==T){## rebuild the heatmap matrix and other datastructures if the flag is set
  cat("Rebuilding matrix\n")
  ## read the gwas results
  cat("Reading the GWAS results...\n")
  d<-read.table(infile, header=T)
  names(d)<-c("chr", "pos","Pvalue", "FRQ","conseq")
  d<-d[!is.na(d$pos),]
  ## align the frequencies to the minor allele
  d$FRQ[d$FRQ > 0.5]<-(1-(d$FRQ[d$FRQ>0.5]))
  
  ## check that the chromosome column is in correct format.
  waitifnot(is.numeric(d$chr), "chr column in gwas data should be numeric, please check if encoded X, or with `chr` prefix")
  ## Support for the Z chromosome?
  lastchr<-max(unique(d$chr))
  
  ## read the snp info
  snp.info<-read.table(snpfile, header=T, sep="\t")
  snp.info<-snp.info[order(snp.info$chr, snp.info$pos, decreasing=F),]
  
  config<-read.table(configfile,sep="\t", header =T,stringsAsFactors = F, skip=10)
  
  pvals<-seq(from=0, to=max.pval, by=pval.split) # max(-log10(d$Pvalue))
  pvals.cells.index<-data.frame(id=1:161,LP=pvals,UP=c(pvals[2:161],max.pval))
  
  final<-matrix(0, nrow = length(pvals), ncol = 0)
  
  chr.matrix.len<-as.data.frame(matrix(nrow=lastchr,ncol=3))
  names(chr.matrix.len)<-c("length", "cumm", "mid")
  
  pos.chunks<-as.data.frame(matrix(nrow=0,ncol=4))
  names(pos.chunks)<-c("posid","chr","s", "e")
  
  pos.idx<-0
  
  idx.count<-vector(mode="numeric",length=length(config$idx)+1)
  max.cellcount<-0
  
  if(debugflag==T){ ## Generate the log file if the flag is set.
    logfile<-file(debugfile,open="a")
    ## columns for the log output table
    cat("chr","st","en","log10p","idx","\n", file=logfile,append = T,sep="\t")
  }
  
  ## locations to be labelled on the plot.
  pos.interest<-data.frame(marker=character,log10pval=numeric,chr=numeric,pos=numeric,col=character)
  
  for(chr in 1:lastchr){
    cat("chromosome",chr,"out of", lastchr,"\r")
    
    chr.slice<-d[d$chr==chr,]
    
    chunks<-seq(from=min(chr.slice$pos), to=max(chr.slice$pos), by=pos.split)
    chunks[length(chunks)+1]<-max(chr.slice$pos)
    
    
    mdat<-matrix(0, nrow = length(pvals), ncol = length(chunks)-1)
    
    for (i in 1:(length(chunks)-1)){
      slice<-chr.slice[chr.slice$pos >= chunks[i] & chr.slice$pos < chunks[i+1],]
      
      for (j in 1:length(pvals)){
        if(j == length(pvals)){ ## last element in pval array - include all variants gt then this
          p.val.slice<-slice[-log10(slice$Pvalue) >= pvals[j],]  
        }
        else{ ## otherwise slice between the current and next values in the pval array
          p.val.slice<-slice[-log10(slice$Pvalue) >= pvals[j] & -log10(slice$Pvalue) < pvals[j+1],]
        }      
        
        len<-dim(p.val.slice)[1]
        idx<-0
        
        conseq.len<-dim(p.val.slice[p.val.slice$conseq==1,])[1] ## test if any high consequences
        maf.len<-dim(p.val.slice[p.val.slice$FRQ<=MAF,])[1] ## test if any MAF < 5%
        
        if(len == 0){idx<-0} ## the case with no variants - blank cell.
        else{
          ## len - number of snps in the bin
          ## conseq.len - does it include any high impact variants
          ## maf.len - does it include any MAF < 5%
          
          ## this defines the region which should be greyed out
          if(pvals[j] <= -log10(FDR)){
            if((chr %% 2) != 0){idx<-config$idx[config$type=="oddchr"]}
            else{idx<-config$idx[config$type=="evenchr"]}
          }
          else{
            for(k in 1:length(config$idx)){
              if(config$type[k] == "val"){ ## Check the config is a valid type.
                  len.chk<-FALSE
                  conseq.chk<-FALSE
                  maf.chk<-FALSE
                  
                  ## check counts
                  if(len >= config$min.count[k]){ ## check that the current cell length is gt or eq than the config min count
                    if(is.na(config$max.count[k])){ ## if the config max is NA then accept length condition
                      len.chk<-TRUE
                    }
                    else if(len < config$max.count[k]){ ## if the current cell length is lt the config max count then accept the length condition
                      len.chk<-TRUE
                    }
                  }
                  
                  ## check HIGH impact
                  if(config$conseq[k] == TRUE && conseq.len > 0){conseq.chk<-TRUE} ## accept if high impact is active in config and there are 1 or more high impact variants within the cell
                  else if(config$conseq[k] == FALSE && conseq.len == 0){conseq.chk<-TRUE} ## accept if high impact is inactive and there are 0 high impact variants within the cell.
                  
                  ## check MAF
                  if(config$maf[k] == TRUE && maf.len > 0){maf.chk<-TRUE} ## accept if MAF 5% active and there is one or more variant with MAF lt 5% in the cell
                  else if(config$maf[k] == FALSE && maf.len == 0){maf.chk<-TRUE} ## accept if MAF 5% active and there are 0 variants with MAF lt 5% in the cell
                  
                  if(len.chk==TRUE && conseq.chk==TRUE && maf.chk==TRUE){ ## if all three clauses are correct then accept the idx for the config and exit the for loop
                    idx<-config$idx[k]
                    
                    if(config$report[k]==TRUE){ ## if reporting is active for the config then add an entry to the pos.interest table.
                      tmp.df<-data.frame(
                        marker=paste(idx,sep=""),
                        log10pval=pvals[j],chr=chr,pos=(chunks[i]+1),
                        col=config$col[k])
                      
                      pos.interest<-rbind(pos.interest,tmp.df)
                    }
                    break ## found config which fulfils cell criteria accept and exit group.
                  }
                  
              }
            }
            
            if(idx == 0){ ## if no idx was found in the annotations, generate a remaining idx (max idx in annotations + 1)
              idx<-max(config$idx)+1; ## assign to others
            }
            
          if(len > max.cellcount){max.cellcount<-len} ## assign the maximum cell count if length of current cell is gt than current
          }
        }

        idx.count[idx]<-(idx.count[idx] + 1) # increment the block count for each index
        
        mdat[j,i]<-idx  
    
        if(debugflag==T){ ## log information for current cell.
          
          oddchridx<-config[config$type=="oddchr",]$idx
          evenchridx<-config[config$type=="evenchr",]$idx
          
          if(idx != oddchridx || idx != evenchridx || idx != 0){ ## generate the log file, dont log the greyed out regions
            cat(chr,chunks[i],chunks[i+1],pvals[j],idx,"\n", file=logfile,append = T,sep="\t")
          }
        }
        
      }
    }
    
    ## cell sizes of the chromosomes
    chr.matrix.len$length[chr]<-dim(mdat)[2]
    chr.matrix.len$cumm[chr]<-dim(mdat)[2]+ dim(final)[2]
    chr.matrix.len$mid[chr]<-chr.matrix.len$cumm[chr] - (chr.matrix.len$length[chr] / 2)
    
    ## bind to the final matrix
    final<-cbind(final,mdat)
    
    tmp.chunks<-as.data.frame(matrix(nrow=length(chunks)-1,ncol=4))
    names(tmp.chunks)<-c("posid","chr","s", "e")
      
    for(j in 1:dim(tmp.chunks)[1]){ ## assign position chunk indexes to table
        #print(j)
        tmp.chunks$posid[j]<-pos.idx+j
        tmp.chunks$chr[j]<-chr
        tmp.chunks$s[j]<-chunks[j]
        tmp.chunks$e[j]<-chunks[j+1]
    }
    
    ## the cell id for each of the position chunks in the heatmap
    pos.chunks<-rbind(pos.chunks,tmp.chunks)
    pos.idx<-pos.idx+dim(tmp.chunks)[1]
  }
  
  if(debugflag==T){ ## close the log file if debuging is active.
    close(logfile) 
  }
  
  snpcells<-vector(length=length(snp.info$markername))
  
  for (i in 1:length(snp.info$markername)){ ## assign a position cell for each of the variants in the SNP information table.
    snpcells[i]<-chrpos.cell(snp.info$chr[i],snp.info$pos[i],pos.chunks)
  }
  
  cat("Melting matrix...\n")
  m<-melt(final)
  names(m)<-c("pval","pos", "val")

  if(dim(pos.interest)[1] > 0){  ## if there are any positions of interest assign the cell positions on the heatmap
    pos.interest$pvalidx<-log10.index(pos.interest$log10pval)
    pos.interest$pos.idx<--1
    
    for (i in 1:length(pos.interest$marker)){ ## assign a position cell for each of the cells of interest
      pos.interest$pos.idx[i]<-chrpos.cell(pos.interest$chr[i],pos.interest$pos[i],pos.chunks)
    }
  }

} ## END OF REBUILD SECTION


## assign the maximum variants in the cells
peak.val<-ceiling(max.cellcount / 100)*100 ## nearest 10,000 to max cell value

col.discrete<-c("white",config$col)

col.text<-vector(mode="character",length=length(config$type))

for(k in 1:length(col.text)){ ## Build the text for the legend using the information in the config table
  if(config$type[k]=="val"){
    
    max.count<-config$max.count[k]
    
    if(is.na(max.count)){ ## if the max count is assigned to NA then the highest cell count (with approp ceiling) is used
      max.count<-peak.val
    }
    else{
      max.count<-(max.count-1)
    }
    
    if(config$min.count[k]==1 &&max.count==1){
      col.text.tmp<-paste(config$idx[k],") ",
                          config$min.count[k],
                          sep="")
    }
    else{
      col.text.tmp<-paste(config$idx[k],") ",
                         config$min.count[k]," - ",
                         max.count,
                         sep="")
    }
  
  if(config$conseq[k]==FALSE && config$maf[k]==FALSE){ ## if both MAF and HIGH impact are disabled - do not add any text
    ##col.text[k]<-col.text.tmp
  }
  else if(config$conseq[k] == TRUE && config$maf[k]==TRUE){ ## if both MAF and HIGH impact are active - print BOTH
    col.text.tmp<-paste(col.text.tmp," (BOTH)",sep="")
  }
  else if(config$conseq[k] == TRUE && config$maf[k] == FALSE){ ## if only HIGH impact is active - print HIGH impact.
    col.text.tmp<-paste(col.text.tmp," (HIGH impact)",sep="")
  }
  else if(config$conseq[k] == FALSE && config$maf[k] == TRUE){ ## if only MAD is active - print MAF
    col.text.tmp<-paste(col.text.tmp," (MAF < ",MAF,")",sep="")
  }
  else{ ## all possible scenarios are covered - should not get here
    cat("Should not get here!\n")  
    stopifnot(FALSE)
  }
    col.text[k]<-paste(col.text.tmp," (",idx.count[k],")",sep="")
  }
  else{ ## otherwise print the config type.
    col.text[k]<-config$type[k]
  }
}

col.brks<-1:length(config$type)

if(max(m$val) > max(config$idx)){ ## if there are any cells which do not have a category add them to the remaining category.
  col.brks<-c(col.brks,max(m$val))
  col.discrete<-c(col.discrete,"orange") ## hard coded colour for remaining.
  col.text<-c(col.text,paste("Remaining (", idx.count[max(m$val)],")",sep=""))
}

pval.seq<-seq(from=pval.units,to=max.pval,by=pval.units)

y.labels<-c("",pval.seq)
y.breaks<-c(0.5,log10.index(pval.seq))

if(showgenes==FALSE){ ## if show genes flag is not set then put all the genes in the table
  snp.info$novel=TRUE
}

snp.info.known<-snp.info[snp.info$novel==FALSE,]
snp.info.novel<-snp.info[snp.info$novel==TRUE,]

main.core<-ggplot(data=m, aes(x=pos,y=pval)) + 
  geom_tile(aes(fill = val))+ ##,colour= val), size=0.01) + 
  theme(legend.position="left",legend.key.size=unit(0.5,"line"),
        legend.title=element_text(size=5),
        legend.text=element_text(size=5)) +
  geom_hline(yintercept=p.val.cell(GWS)+0.5, linetype="dashed") + ## GWS line
  geom_hline(yintercept=p.val.cell(FDR)+0.5, linetype="dashed") + ## FDR line
  scale_fill_gradientn(colours = col.discrete, 
        guide="legend", breaks=col.brks, 
        labels=col.text,name = "Variant Count") + 
  scale_y_continuous("-log10(p)",labels=y.labels, 
        breaks=y.breaks,
        expand=c(0,0),trans="reverse",position="right") +
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank())   + coord_flip(expand=T) + 
  theme(plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_x_continuous("", labels=rep("",(lastchr*2)), 
        breaks=c(chr.matrix.len$mid,chr.matrix.len$cumm),
        expand=c(0,0),trans="reverse",position="top") +
  theme(axis.ticks.y = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5))
  
if(dim(snp.info.known)[1] > 0){ ## if there is one or more known SNPs in the table then label the manhattan plot with them.
  repel.df<-as.data.frame(matrix(nrow=dim(snp.info.known)[1],ncol=3))
  names(repel.df)<-c("marker","pvalidx","posidx")
  
  repel.df$marker<-snp.info.known$NearestGene
  if(showrsids==T){
    repel.df$marker<-snp.info.known$markername
  }
  
  repel.df$pvalidx<-p.val.index(snp.info.known$Pvalue) 
  repel.df$posidx<-snpcells[snp.info$novel==FALSE] 
  
  repel.df[repel.df$pvalidx > log10.index(max.pval),]$pvalidx<-log10.index(max.pval)
  
  final.repel.plot<-main.core+
    geom_label_repel(data=repel.df, aes(posidx,pvalidx, label=marker),
          size=textsize, force=1, nudge_y = 10,nudge_x=10,
          segment.colour="black", min.segment.length = 0,
          segment.size=0.25, seed=500, max.iter = 5000,
          point.padding = NA)
}

if(dim(pos.interest)[1]> 0){ ## if any positions of interest label on the manhattan plot then label them
  main.core<-main.core+ geom_label_repel(data=pos.interest, aes(pos.idx,pvalidx, label=marker),
        size=textsize, force=5, nudge_y = 10,nudge_x=10,
        segment.colour="black", min.segment.length = 0,
        segment.size=0.25, seed=500, max.iter = 5000,
        point.padding = NA,segment.color = "black",color="black")
}

title.pos<-c(log10.index(17)+1, log10.index(13)+1, log10.index(11.5)+1,log10.index(10)+1,log10.index(7)+1)

## Generate table for the novel genes
## Start with a completely blank plot (table1)
table1<-ggplot(data=m, aes(x=pos,y=pval)) + 
  geom_tile(aes(fill = rep(0,dim(m)[1]))) +
  scale_x_continuous("", labels=rep("",(lastchr*2)),
        breaks=c(chr.matrix.len$mid,chr.matrix.len$cumm),
        expand=c(0,0),trans="reverse",position="top") + #, minor_breaks=chr.matrix.len$cumm) +
  scale_y_continuous("",labels=rep("",length(y.labels)),
        breaks=y.breaks,
        expand=c(0,0),trans="reverse",position="right") +
  coord_flip(expand=T) + 
  scale_fill_gradientn(colours = c("white","white"), 
        guide=FALSE, breaks=c(0,1), 
        labels=c("0","1"),name = "")  +
  theme(axis.ticks.y = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))


lim<-layer_scales(table1)
xmin<-lim$x$range$range[2]*-1
xmax<-lim$x$range$range[1]*-1

segment.indexes<-c(log10.index(17.5)+1,log10.index(19)+1,log10.index(20)+1)
chr.num.pos<-(log10.index(19.75)+1)
brks.pos<-(log10.index(19.7)+1)

text.pos<-seq(from=xmin, to=xmax,length.out=length(snp.info.novel$markername)+1)
text.pos1<-text.pos[2:length(text.pos)]
title.pos1<-text.pos[1]

## add the SNP information to the table
## Can add additional columns to the table here.
table2<-table1+
  annotate("text", x = text.pos1,
           y = title.pos[1], label = as.character(snp.info.novel$markername),
           angle=0,size=textsize, hjust=0) + ## markername
  annotate("text", x = text.pos1, 
           y = title.pos[2], label = as.character(format(round(snp.info.novel$eaf,2),nsmall=2)),
           angle=0,size=textsize, hjust=0) + ## EAF
  annotate("text", x = text.pos1, 
           y = title.pos[3], label = as.character(format(round(snp.info.novel$OR,2),nsmall=2)),
           angle=0,size=textsize, hjust=0) + ## OR
  annotate("text", x = text.pos1, 
           y = title.pos[4], label = as.character(formatC(snp.info.novel$Pvalue, format = "E", digits = 2)),
           angle=0,size=textsize, hjust=0) + ## P-value
  annotate("text", x = text.pos1,
           y = title.pos[5], label =as.character(snp.info.novel$NearestGene),
           angle=0,size=textsize, hjust=0,fontface = 'italic') + ## Nearest Gene
  annotate("segment", x = text.pos1,
           xend = snpcells[snp.info$novel==TRUE], y = segment.indexes[1], yend = segment.indexes[2],
           colour = "blue", linetype="dashed", size=0.5) +  ## segment from midpoint to table row
  annotate("segment", x = snpcells[snp.info$novel==TRUE], 
           xend = snpcells[snp.info$novel==TRUE], y = segment.indexes[2], yend = segment.indexes[3],
           colour = "blue", linetype="dashed", size=0.5) + ## segment from axis to mid point
  annotate("text", x = chr.matrix.len$mid,
           y = chr.num.pos, label = as.character(1:lastchr), 
           angle=0,size=3.5, hjust=0)+ ## chromosome labels
  annotate("segment", x = chr.matrix.len$cumm, 
           xend = chr.matrix.len$cumm, y = brks.pos, yend = segment.indexes[3], 
           colour = "black", linetype="solid") + ## x axis breaks
  annotate("segment", x = 0,  
           xend = xmax, y = segment.indexes[3], yend = segment.indexes[3],
           colour = "black", linetype="solid")  ## xaxis solid line

## add the title to the table.
## Can add additional table headers here.
final.table.plot<-table2 + 
  annotate("text", x = title.pos1,
           y = title.pos[1], label = "SNP",
           angle=0,size=textsize, hjust=0,fontface = 'bold') +
  annotate("text", x = title.pos1,
           y = title.pos[2], label = "EAF",
           angle=0,size=textsize, hjust=0,fontface = 'bold') +
  annotate("text", x = title.pos1,
           y = title.pos[3], label = "OR",
           angle=0,size=textsize, hjust=0,fontface = 'bold') +
  annotate("text", x = title.pos1,
           y = title.pos[4], label = "p-value",
           angle=0,size=textsize, hjust=0,fontface = 'bold') +
  annotate("text", x = title.pos1,
           y = title.pos[5], label = "Gene",
           angle=0,size=textsize, hjust=0,fontface = 'bold') 

gt <- ggplot_gtable(ggplot_build(final.table.plot)) # p4
gt$layout$clip[gt$layout$name == "panel"] <- "off"

## hard code variables for positions of two plots on qplot
manh.max<-7
annot.min<-6.7

## always draw as a PDF.

if(drawastiff==T){
  tiff(filename = paste(outfile,".tif",sep=""),width = 8.27,height = 11.69, units="in",res=300)
} else{
  pdf(paste(outfile,".pdf",sep=""),width = 8.27,height = 11.69,onefile = F)
}




final.plot<-main.core
 
if(showgenes==TRUE){ ## if show genes flag is set then assign to the final plot.
  final.plot<-final.repel.plot
}

qplot(1:10,1:10,colour=I("white")) +
  annotation_custom(grob=ggplotGrob(final.plot), xmin=0.5,xmax=manh.max, ymin=1,ymax = 10) +
  annotation_custom(grob=gt, xmin=annot.min, xmax=10.5, ymin=1,ymax = 10) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(plot.background = element_rect(fill = 'white', colour = 'white')) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))

dev.off()

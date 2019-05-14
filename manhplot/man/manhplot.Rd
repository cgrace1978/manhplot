\name{manhplot}
\alias{manhplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
MANHATTAN++
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
manhplot(infile, outfile, configfile, snpfile, drawastiff = F, GWS = 5e-08, FDR = 0.001, MAF = 0.05, chrname = "chr", posname = "pos", pvalname = "pvalue", frqname = "maf", conseqname = "conseq")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{infile}{
%%     ~~Describe \code{infile} here~~
}
  \item{outfile}{
%%     ~~Describe \code{outfile} here~~
}
  \item{configfile}{
%%     ~~Describe \code{configfile} here~~
}
  \item{snpfile}{
%%     ~~Describe \code{snpfile} here~~
}
  \item{drawastiff}{
%%     ~~Describe \code{drawastiff} here~~
}
  \item{GWS}{
%%     ~~Describe \code{GWS} here~~
}
  \item{FDR}{
%%     ~~Describe \code{FDR} here~~
}
  \item{MAF}{
%%     ~~Describe \code{MAF} here~~
}
  \item{chrname}{
%%     ~~Describe \code{chrname} here~~
}
  \item{posname}{
%%     ~~Describe \code{posname} here~~
}
  \item{pvalname}{
%%     ~~Describe \code{pvalname} here~~
}
  \item{frqname}{
%%     ~~Describe \code{frqname} here~~
}
  \item{conseqname}{
%%     ~~Describe \code{conseqname} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (infile, outfile, configfile, snpfile, drawastiff = F, 
    GWS = 5e-08, FDR = 0.001, MAF = 0.05, chrname = "chr", posname = "pos", 
    pvalname = "pvalue", frqname = "maf", conseqname = "conseq") 
{
    requireNamespace("reshape2", quietly = TRUE)
    requireNamespace("ggplot2", quietly = TRUE)
    requireNamespace("ggrepel", quietly = TRUE)
    requireNamespace("gridExtra", quietly = TRUE)
    showgenes <- F
    showrsids <- F
    pos.split <- 3e+06
    pval.split <- 0.125
    max.pval <- 20
    pval.units <- 5
    textsize <- 2
    rebuild <- T
    debugfile <- paste("manh_", format(Sys.time(), "\%d-\%m-\%y.\%H-\%M-\%S"), 
        ".log", sep = "")
    debugflag <- F
    waitifnot <- function(cond, mess) {
        if (!cond) {
            message(mess)
            message(deparse(substitute(cond)), " is not TRUE")
            while (TRUE) {
            }
        }
    }
    log10.index <- function(log) {
        idx <- (log/pval.split) + 0.5
        return(idx)
    }
    p.val.index <- function(p.val) {
        gws <- -1 * log10(p.val)
        idx <- log10.index(gws)
        return(idx)
    }
    p.val.cell <- function(p.val, pval.chunks) {
        gws <- -1 * log10(p.val)
        idx <- log10.cell(gws)
        return(idx)
    }
    log10.cell <- function(log, pval.chunks) {
        idx <- pvals.cells.index$id[log >= pvals.cells.index$LP & 
            log < pvals.cells.index$UP]
        return(idx)
    }
    chrpos.cell <- function(chr, position, chr.chunks) {
        slice <- chr.chunks[chr.chunks$chr == chr & chr.chunks$s < 
            position & chr.chunks$e > position, ]
        return(slice$posid)
    }
    if (rebuild == T) {
        cat("Rebuilding matrix\n")
        cat("Reading the GWAS results...\n")
        d <- read.table(infile, header = T)
        names <- names(d)
        for (i in 1:length(names)) {
            col <- names[i]
            if (col == chrname) {
                names[i] <- "chr"
            }
            if (col == posname) {
                names[i] <- "pos"
            }
            if (col == pvalname) {
                names[i] <- "Pvalue"
            }
            if (col == frqname) {
                names[i] <- "FRQ"
            }
            if (col == conseqname) {
                names[i] <- "conseq"
            }
        }
        names(d) <- names
        d <- d[!is.na(d$pos), ]
        d$FRQ[d$FRQ > 0.5] <- (1 - (d$FRQ[d$FRQ > 0.5]))
        waitifnot(is.numeric(d$chr), "chr column in gwas data should be numeric, please check if encoded X, or with `chr` prefix")
        lastchr <- max(unique(d$chr))
        snp.info <- read.table(snpfile, header = T, sep = "\t")
        snp.info <- snp.info[order(snp.info$chr, snp.info$pos, 
            decreasing = F), ]
        config <- read.table(configfile, sep = "\t", header = T, 
            stringsAsFactors = F, skip = 10)
        pvals <- seq(from = 0, to = max.pval, by = pval.split)
        pvals.cells.index <- data.frame(id = 1:161, LP = pvals, 
            UP = c(pvals[2:161], max.pval))
        final <- matrix(0, nrow = length(pvals), ncol = 0)
        chr.matrix.len <- as.data.frame(matrix(nrow = lastchr, 
            ncol = 3))
        names(chr.matrix.len) <- c("length", "cumm", "mid")
        pos.chunks <- as.data.frame(matrix(nrow = 0, ncol = 4))
        names(pos.chunks) <- c("posid", "chr", "s", "e")
        pos.idx <- 0
        idx.count <- vector(mode = "numeric", length = length(config$idx) + 
            1)
        max.cellcount <- 0
        if (debugflag == T) {
            logfile <- file(debugfile, open = "a")
            cat("chr", "st", "en", "log10p", "idx", "\n", file = logfile, 
                append = T, sep = "\t")
        }
        pos.interest <- data.frame(marker = character, log10pval = numeric, 
            chr = numeric, pos = numeric, col = character)
        for (chr in 1:lastchr) {
            cat("chromosome", chr, "out of", lastchr, "\r")
            chr.slice <- d[d$chr == chr, ]
            chunks <- seq(from = min(chr.slice$pos), to = max(chr.slice$pos), 
                by = pos.split)
            chunks[length(chunks) + 1] <- max(chr.slice$pos)
            mdat <- matrix(0, nrow = length(pvals), ncol = length(chunks) - 
                1)
            for (i in 1:(length(chunks) - 1)) {
                slice <- chr.slice[chr.slice$pos >= chunks[i] & 
                  chr.slice$pos < chunks[i + 1], ]
                for (j in 1:length(pvals)) {
                  if (j == length(pvals)) {
                    p.val.slice <- slice[-log10(slice$Pvalue) >= 
                      pvals[j], ]
                  }
                  else {
                    p.val.slice <- slice[-log10(slice$Pvalue) >= 
                      pvals[j] & -log10(slice$Pvalue) < pvals[j + 
                      1], ]
                  }
                  len <- dim(p.val.slice)[1]
                  idx <- 0
                  conseq.len <- dim(p.val.slice[p.val.slice$conseq == 
                    1, ])[1]
                  maf.len <- dim(p.val.slice[p.val.slice$FRQ <= 
                    MAF, ])[1]
                  if (len == 0) {
                    idx <- 0
                  }
                  else {
                    if (pvals[j] <= -log10(FDR)) {
                      if ((chr\%\%2) != 0) {
                        idx <- config$idx[config$type == "oddchr"]
                      }
                      else {
                        idx <- config$idx[config$type == "evenchr"]
                      }
                    }
                    else {
                      for (k in 1:length(config$idx)) {
                        if (config$type[k] == "val") {
                          len.chk <- FALSE
                          conseq.chk <- FALSE
                          maf.chk <- FALSE
                          if (len >= config$min.count[k]) {
                            if (is.na(config$max.count[k])) {
                              len.chk <- TRUE
                            }
                            else if (len < config$max.count[k]) {
                              len.chk <- TRUE
                            }
                          }
                          if (config$conseq[k] == TRUE && conseq.len > 
                            0) {
                            conseq.chk <- TRUE
                          }
                          else if (config$conseq[k] == FALSE && 
                            conseq.len == 0) {
                            conseq.chk <- TRUE
                          }
                          if (config$maf[k] == TRUE && maf.len > 
                            0) {
                            maf.chk <- TRUE
                          }
                          else if (config$maf[k] == FALSE && 
                            maf.len == 0) {
                            maf.chk <- TRUE
                          }
                          if (len.chk == TRUE && conseq.chk == 
                            TRUE && maf.chk == TRUE) {
                            idx <- config$idx[k]
                            if (config$report[k] == TRUE) {
                              tmp.df <- data.frame(marker = paste(idx, 
                                sep = ""), log10pval = pvals[j], 
                                chr = chr, pos = (chunks[i] + 
                                  1), col = config$col[k])
                              pos.interest <- rbind(pos.interest, 
                                tmp.df)
                            }
                            break
                          }
                        }
                      }
                      if (idx == 0) {
                        idx <- max(config$idx) + 1
                      }
                      if (len > max.cellcount) {
                        max.cellcount <- len
                      }
                    }
                  }
                  idx.count[idx] <- (idx.count[idx] + 1)
                  mdat[j, i] <- idx
                  if (debugflag == T) {
                    oddchridx <- config[config$type == "oddchr", 
                      ]$idx
                    evenchridx <- config[config$type == "evenchr", 
                      ]$idx
                    if (idx != oddchridx || idx != evenchridx || 
                      idx != 0) {
                      cat(chr, chunks[i], chunks[i + 1], pvals[j], 
                        idx, "\n", file = logfile, append = T, 
                        sep = "\t")
                    }
                  }
                }
            }
            chr.matrix.len$length[chr] <- dim(mdat)[2]
            chr.matrix.len$cumm[chr] <- dim(mdat)[2] + dim(final)[2]
            chr.matrix.len$mid[chr] <- chr.matrix.len$cumm[chr] - 
                (chr.matrix.len$length[chr]/2)
            final <- cbind(final, mdat)
            tmp.chunks <- as.data.frame(matrix(nrow = length(chunks) - 
                1, ncol = 4))
            names(tmp.chunks) <- c("posid", "chr", "s", "e")
            for (j in 1:dim(tmp.chunks)[1]) {
                tmp.chunks$posid[j] <- pos.idx + j
                tmp.chunks$chr[j] <- chr
                tmp.chunks$s[j] <- chunks[j]
                tmp.chunks$e[j] <- chunks[j + 1]
            }
            pos.chunks <- rbind(pos.chunks, tmp.chunks)
            pos.idx <- pos.idx + dim(tmp.chunks)[1]
        }
        if (debugflag == T) {
            close(logfile)
        }
        snpcells <- vector(length = length(snp.info$markername))
        for (i in 1:length(snp.info$markername)) {
            snpcells[i] <- chrpos.cell(snp.info$chr[i], snp.info$pos[i], 
                pos.chunks)
        }
        cat("\nMelting matrix...\n")
        m <- melt(final)
        names(m) <- c("pval", "pos", "val")
        if (dim(pos.interest)[1] > 0) {
            pos.interest$pvalidx <- log10.index(pos.interest$log10pval)
            pos.interest$pos.idx <- -1
            for (i in 1:length(pos.interest$marker)) {
                pos.interest$pos.idx[i] <- chrpos.cell(pos.interest$chr[i], 
                  pos.interest$pos[i], pos.chunks)
            }
        }
    }
    peak.val <- ceiling(max.cellcount/100) * 100
    col.discrete <- c("white", config$col)
    col.text <- vector(mode = "character", length = length(config$type))
    for (k in 1:length(col.text)) {
        if (config$type[k] == "val") {
            max.count <- config$max.count[k]
            if (is.na(max.count)) {
                max.count <- peak.val
            }
            else {
                max.count <- (max.count - 1)
            }
            if (config$min.count[k] == 1 && max.count == 1) {
                col.text.tmp <- paste(config$idx[k], ") ", config$min.count[k], 
                  sep = "")
            }
            else {
                col.text.tmp <- paste(config$idx[k], ") ", config$min.count[k], 
                  " - ", max.count, sep = "")
            }
            if (config$conseq[k] == FALSE && config$maf[k] == 
                FALSE) {
            }
            else if (config$conseq[k] == TRUE && config$maf[k] == 
                TRUE) {
                col.text.tmp <- paste(col.text.tmp, " (BOTH)", 
                  sep = "")
            }
            else if (config$conseq[k] == TRUE && config$maf[k] == 
                FALSE) {
                col.text.tmp <- paste(col.text.tmp, " (HIGH impact)", 
                  sep = "")
            }
            else if (config$conseq[k] == FALSE && config$maf[k] == 
                TRUE) {
                col.text.tmp <- paste(col.text.tmp, " (MAF < ", 
                  MAF, ")", sep = "")
            }
            else {
                cat("Should not get here!\n")
                stopifnot(FALSE)
            }
            col.text[k] <- paste(col.text.tmp, " (", idx.count[k], 
                ")", sep = "")
        }
        else {
            col.text[k] <- config$type[k]
        }
    }
    col.brks <- 1:length(config$type)
    if (max(m$val) > max(config$idx)) {
        col.brks <- c(col.brks, max(m$val))
        col.discrete <- c(col.discrete, "orange")
        col.text <- c(col.text, paste("Remaining (", idx.count[max(m$val)], 
            ")", sep = ""))
    }
    pval.seq <- seq(from = pval.units, to = max.pval, by = pval.units)
    y.labels <- c("", pval.seq)
    y.breaks <- c(0.5, log10.index(pval.seq))
    if (showgenes == FALSE) {
        snp.info$novel = TRUE
    }
    snp.info.known <- snp.info[snp.info$novel == FALSE, ]
    snp.info.novel <- snp.info[snp.info$novel == TRUE, ]
    main.core <- ggplot(data = m, aes(x = pos, y = pval)) + geom_tile(aes(fill = val)) + 
        theme(legend.position = "left", legend.key.size = unit(0.5, 
            "line"), legend.title = element_text(size = 5), legend.text = element_text(size = 5)) + 
        geom_hline(yintercept = p.val.cell(GWS) + 0.5, linetype = "dashed") + 
        geom_hline(yintercept = p.val.cell(FDR) + 0.5, linetype = "dashed") + 
        scale_fill_gradientn(colours = col.discrete, guide = "legend", 
            breaks = col.brks, labels = col.text, name = "Variant Count") + 
        scale_y_continuous("-log10(p)", labels = y.labels, breaks = y.breaks, 
            expand = c(0, 0), trans = "reverse", position = "right") + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        coord_flip(expand = T) + theme(plot.margin = unit(c(0, 
        0, 0, 0), "cm")) + scale_x_continuous("", labels = rep("", 
        (lastchr * 2)), breaks = c(chr.matrix.len$mid, chr.matrix.len$cumm), 
        expand = c(0, 0), trans = "reverse", position = "top") + 
        theme(axis.ticks.y = element_blank()) + theme(axis.line.x = element_line(color = "black", 
        size = 0.5))
    if (dim(snp.info.known)[1] > 0) {
        repel.df <- as.data.frame(matrix(nrow = dim(snp.info.known)[1], 
            ncol = 3))
        names(repel.df) <- c("marker", "pvalidx", "posidx")
        repel.df$marker <- snp.info.known$NearestGene
        if (showrsids == T) {
            repel.df$marker <- snp.info.known$markername
        }
        repel.df$pvalidx <- p.val.index(snp.info.known$Pvalue)
        repel.df$posidx <- snpcells[snp.info$novel == FALSE]
        repel.df[repel.df$pvalidx > log10.index(max.pval), ]$pvalidx <- log10.index(max.pval)
        final.repel.plot <- main.core + geom_label_repel(data = repel.df, 
            aes(posidx, pvalidx, label = marker), size = textsize, 
            force = 1, nudge_y = 10, nudge_x = 10, segment.colour = "black", 
            min.segment.length = 0, segment.size = 0.25, seed = 500, 
            max.iter = 5000, point.padding = NA)
    }
    if (dim(pos.interest)[1] > 0) {
        main.core <- main.core + geom_label_repel(data = pos.interest, 
            aes(pos.idx, pvalidx, label = marker), size = textsize, 
            force = 5, nudge_y = 10, nudge_x = 10, segment.colour = "black", 
            min.segment.length = 0, segment.size = 0.25, seed = 500, 
            max.iter = 5000, point.padding = NA, segment.color = "black", 
            color = "black")
    }
    title.pos <- c(log10.index(17) + 1, log10.index(13) + 1, 
        log10.index(11.5) + 1, log10.index(10) + 1, log10.index(7) + 
            1)
    table1 <- ggplot(data = m, aes(x = pos, y = pval)) + geom_tile(aes(fill = rep(0, 
        dim(m)[1]))) + scale_x_continuous("", labels = rep("", 
        (lastchr * 2)), breaks = c(chr.matrix.len$mid, chr.matrix.len$cumm), 
        expand = c(0, 0), trans = "reverse", position = "top") + 
        scale_y_continuous("", labels = rep("", length(y.labels)), 
            breaks = y.breaks, expand = c(0, 0), trans = "reverse", 
            position = "right") + coord_flip(expand = T) + scale_fill_gradientn(colours = c("white", 
        "white"), guide = FALSE, breaks = c(0, 1), labels = c("0", 
        "1"), name = "") + theme(axis.ticks.y = element_blank()) + 
        theme(axis.ticks.x = element_blank()) + theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.border = element_blank(), 
        panel.background = element_blank()) + theme(plot.margin = unit(c(0, 
        0, 0, 0), "cm"))
    lim <- layer_scales(table1)
    xmin <- lim$x$range$range[2] * -1
    xmax <- lim$x$range$range[1] * -1
    segment.indexes <- c(log10.index(17.5) + 1, log10.index(19) + 
        1, log10.index(20) + 1)
    chr.num.pos <- (log10.index(19.75) + 1)
    brks.pos <- (log10.index(19.7) + 1)
    text.pos <- seq(from = xmin, to = xmax, length.out = length(snp.info.novel$markername) + 
        1)
    text.pos1 <- text.pos[2:length(text.pos)]
    title.pos1 <- text.pos[1]
    table2 <- table1 + annotate("text", x = text.pos1, y = title.pos[1], 
        label = as.character(snp.info.novel$markername), angle = 0, 
        size = textsize, hjust = 0) + annotate("text", x = text.pos1, 
        y = title.pos[2], label = as.character(format(round(snp.info.novel$eaf, 
            2), nsmall = 2)), angle = 0, size = textsize, hjust = 0) + 
        annotate("text", x = text.pos1, y = title.pos[3], label = as.character(format(round(snp.info.novel$OR, 
            2), nsmall = 2)), angle = 0, size = textsize, hjust = 0) + 
        annotate("text", x = text.pos1, y = title.pos[4], label = as.character(formatC(snp.info.novel$Pvalue, 
            format = "E", digits = 2)), angle = 0, size = textsize, 
            hjust = 0) + annotate("text", x = text.pos1, y = title.pos[5], 
        label = as.character(snp.info.novel$NearestGene), angle = 0, 
        size = textsize, hjust = 0, fontface = "italic") + annotate("segment", 
        x = text.pos1, xend = snpcells[snp.info$novel == TRUE], 
        y = segment.indexes[1], yend = segment.indexes[2], colour = "blue", 
        linetype = "dashed", size = 0.5) + annotate("segment", 
        x = snpcells[snp.info$novel == TRUE], xend = snpcells[snp.info$novel == 
            TRUE], y = segment.indexes[2], yend = segment.indexes[3], 
        colour = "blue", linetype = "dashed", size = 0.5) + annotate("text", 
        x = chr.matrix.len$mid, y = chr.num.pos, label = as.character(1:lastchr), 
        angle = 0, size = 3.5, hjust = 0) + annotate("segment", 
        x = chr.matrix.len$cumm, xend = chr.matrix.len$cumm, 
        y = brks.pos, yend = segment.indexes[3], colour = "black", 
        linetype = "solid") + annotate("segment", x = 0, xend = xmax, 
        y = segment.indexes[3], yend = segment.indexes[3], colour = "black", 
        linetype = "solid")
    final.table.plot <- table2 + annotate("text", x = title.pos1, 
        y = title.pos[1], label = "SNP", angle = 0, size = textsize, 
        hjust = 0, fontface = "bold") + annotate("text", x = title.pos1, 
        y = title.pos[2], label = "EAF", angle = 0, size = textsize, 
        hjust = 0, fontface = "bold") + annotate("text", x = title.pos1, 
        y = title.pos[3], label = "OR", angle = 0, size = textsize, 
        hjust = 0, fontface = "bold") + annotate("text", x = title.pos1, 
        y = title.pos[4], label = "p-value", angle = 0, size = textsize, 
        hjust = 0, fontface = "bold") + annotate("text", x = title.pos1, 
        y = title.pos[5], label = "Gene", angle = 0, size = textsize, 
        hjust = 0, fontface = "bold")
    gt <- ggplot_gtable(ggplot_build(final.table.plot))
    gt$layout$clip[gt$layout$name == "panel"] <- "off"
    manh.max <- 7
    annot.min <- 6.7
    if (drawastiff == T) {
        tiff(filename = paste(outfile, ".tif", sep = ""), width = 8.27, 
            height = 11.69, units = "in", res = 300)
    }
    else {
        pdf(paste(outfile, ".pdf", sep = ""), width = 8.27, height = 11.69, 
            onefile = F)
    }
    final.plot <- main.core
    if (showgenes == TRUE) {
        final.plot <- final.repel.plot
    }
    print(qplot(1:10, 1:10, colour = I("white")) + annotation_custom(grob = ggplotGrob(final.plot), 
        xmin = 0.5, xmax = manh.max, ymin = 1, ymax = 10) + annotation_custom(grob = gt, 
        xmin = annot.min, xmax = 10.5, ymin = 1, ymax = 10) + 
        theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
            axis.ticks.x = element_blank()) + theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        theme(plot.background = element_rect(fill = "white", 
            colour = "white")) + theme(panel.background = element_rect(fill = "white", 
        colour = "white")))
    dev.off()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

# QDNAseq.chm13v2: QDNAseq bin annotations for T2T CHM13v2

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R build
status](https://github.com/r-lib/usethis/workflows/R-CMD-check/badge.svg)](https://github.com/r-lib/usethis/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/QDNAseq.chm13v2)](https://CRAN.R-project.org/package=QDNAseq.chm13v2)
<!-- badges: end -->

We provide bin indices to use with QDNAseq at 10, 20, 50, 100, 150, 200, and 500 Kb for the T2T CHM13v2 assembly.  The annotations were created based on the steps from the [QDNAseq vignette](https://bioconductor.org/packages/release/bioc/html/QDNAseq.html), and [QDNAseq.hg38](https://github.com/asntech/QDNAseq.hg38).

**NOTE:  Bins are pending variance estimation with 1000 Genomes Samples. Genomes are currently being downloaded and aligned**

## Installation

You can install the development version of QDNAseq.chm13v2 like so:

``` r
# devtools::install_github("QDNAseq.chm13v2", dependencies = TRUE)
```

## Steps to generate bins

### 1. Prepare T2T CHM13v2

Download T2T assembly [hs1](https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/) and the [GRCh38.d1.vd1.fa](https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files) from TCGA-GDC assembly to extract viral genomes.

```bash
mkdir -p T2T.chm13v2/Sequence/WholeGenomeFasta/  T2T.chm13v2/Annotation/Mappability/
cd T2T.chm13v2/Sequence/WholeGenomeFasta/

## download t2t directory
rsync -avz rsync://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/ .

zcat hs1.fa.gz | sed -e 's/^chr//' | tr '[actg] [ACTG]' > t2t_chm13v2.vd1.fa
cat gdc_viral.fa >> t2t_chm13v2.vd1.fa

```

### 2. Create mappability files
With advances in sequencing, and the growing number of archives we prepared mapabilities at 3 read lengths, 50, 100, and 150 bp using GenMap.  Mappability using GEM was also computed, however, these were highly correlated, thus opted for GenMap as it appears to be maintained.

```bash
cd ../../Annotation/Mappability/
ln -s ../../Sequence/WholeGenomeFasta/t2t_chm13v2.vd1.fa

genmap -F chm13v2.vd1.fa -I chm13v2.vd1.genmap.index

## estimate mappability for 50 100 and 150 bp
for i in 50 100 150 ; do
    [[ !- t2t_genmap_k${i}_E2/ ]] || mkdir -p t2t_genmap_k${i}_E2/
    ##
    genmap -T 12 -E2 -k $i --wig \
           --index chm13v2.vd1.genmap.index/  \
           --output t2t_genmap_k${i}_E2/chm13v2.vd1.k${i}_E2.genmap
done
```

### 3. Excluded Regions
Exclude regions were obtained from R/Bioconductor [ExcludeRanges](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10126321/) package, [here](https://drive.google.com/drive/folders/1sF9m8Y3eZouTZ3IEEywjs2kfHOWFBSJT).

## 4. Preparing Bins
Using [QDNAseq.hg38](https://github.com/asntech/QDNAseq.hg38/tree/main) as a template for binc reation

``` r
## code to generate bins on T2T CHM13v2
## use R/R-4.1.2
## module load R/R-4.1.2

## libraries
library(Biobase)
library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)
library(QDNAseq)
library(future)
library(tidyverse)

#set virtual mem
options(future.globals.maxSize= 7e+10) #8912896000)

# change the current plan to access parallelization
future::plan("multicore", workers = 6)

## ---- prepare bins at various sizes

bin.sizes <- c(1000, 500, 200, 150, 100, 50, 30, 20, 15, 10, 5, 1)
names(bin.sizes) = paste0(bin.sizes, "kbp")

kmers <- c(50, 100, 150)
names(kmers) <- paste0("PE", kmers)
names(kmers)[1] <- "SR50"

t2t.mappability <- "../Mappability/"


## ---- loop through SR50, PE100 and PE150 
for(k in names(kmers)) {

    ## loop through bin.sizes
    for(binsize in bin.sizes) {

        ## create bins
        bins <- createBins(bsgenome = BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0,
                           binSize = binsize, ignoreMitochondria = FALSE)

        ## write out bed file
        bed.file <- paste0("chm13v2.", binsize, "kbp.bed")
        write.table(cbind(bins[,1:3], rownames(bins)),
                    file = bed.file,
                    sep = "\t",
                    row.names = FALSE, col.names = FALSE, quote = FALSE)

        ## mappability files
        t2t.mp.bigwig <- file.path(
            t2t.mappability,
            paste0("t2t_genmap_k", kmers[k],
                   "_E2/chm13v2.vd1.k", kmers[k], "_E2.genmap.bw"))
        
			## exclude list file
        t2t.exclude.file <- file.path(t2t.mappability,
                                      "T2T.excluderanges.bed")
        bigWigAvgOB <- file.path("bigWigAverageOverBed")

    
        ## calculateMappability was not reading k100 nor k150, thus
        ## computing overlaps manually, and importing directly

        cmd <- paste(bigWigAvgOB, t2t.mp.bigwig, bed.file, "stdout | cut -f 1,5 ")

        mappability <- system(cmd, intern = TRUE)
        mappability <- as.data.frame(
            do.call(rbind,
                    strsplit(mappability, "\t"))) %>%
            transform(row.names = V1, V1 = NULL) %>%
            rename("mappability" = "V2")
	bins$mappability <- mappability[rownames(bins), "mappability"]

        ## estimate % overap to exclude list
        bins$blacklist <- calculateBlacklist(
            bins,
            bedFiles = t2t.exclude.file)

        ## make empty residual column
        bins$residual <- NA

        ## make use column
        bins$use <- bins$chromosome %in% as.character(1:22) & bins$bases > 0
		
		## PENDING 
        ## Count bins across 1000 Genomes bam files
        ## tg <- binReadCounts(bins,
		##	                   path = file.path("../1000Genomes/T2T_BIN_VARIANCE/bam_out", k), 
		##                     cache=TRUE)
        ##
        ## bins$residual <- iterateResiduals(tg)

        bins <- AnnotatedDataFrame(
            bins,
            varMetadata = data.frame(
                labelDescription=c(
                    "Chromosome name",
                    "Base pair start position",
                    "Base pair end position",
                    "Percentage of non-N nucleotides (of full bin size)",
                    "Percentage of C and G nucleotides (of non-N nucleotides)",
                    paste0("Average mappability of ", kmers[k],
                           "mers with a maximum of 2 mismatches"),
                    "Percent overlap with ExcludeRanges T2T excluded regions",
                    "Median loess residual from 1000 Genomes (50mers)",
                    "Whether the bin should be used in subsequent analysis steps"),
                row.names = colnames(bins)))

        QDNAseqInfo <- list(
            author="Rodrigo Gularte Merida",
            date=Sys.time(),
            organism='Hsapiens',
            build='CHM13v2',
            version=packageVersion("QDNAseq"),
            url=paste0(
                "https://github.com/RodrigoGM/QDNAseq.chm13v2/raw/master/data/chm13v2.",
                binsize, "kbp.", k, ".rda"),
            md5=digest::digest(bins@data),
            sessionInfo=sessionInfo())

        attr(bins, "QDNAseq") <- QDNAseqInfo
        save(bins, file=file.path("data",
                                  paste0("chm13v2.", binsize, "kbp.", k,".rda")),
             compress='xz')

    }
}

```

The following bash routines were used to change object names from `bins` to match the data, and to document each data set

```bash
#!/bin/bash

bin_name_change() {
	echo "## $1"
	echo "data($1)"
	echo "$1 <- bins"
	echo "usethis::use_data($1, overwrite = TRUE, compress = \"xz\")"
	echo "rm(bins)"
	echo ""
}

for i in $(ls data/ | sed -e 's/.rda//' | sort -k1V) ; do 
	bin_name_change $i 
done > R/name_change.R

R --vanilla -q < R/name_change.R


bin_data_documentation() {
	echo "#' $1"
	echo "#'"
	echo "#' Bin annotations for ${KBP}kbp using $KMER reads"
	echo "#'"
	echo "#' @docType data"
	echo "#' @usage data($1)"
		echo "#' @format A annotated data frame with pre-calcuated metrics for each bin describing:"
	echo "#' \describe{"
	echo "#'   \item{chromosome}{Chromosome name}"
	echo "#'   \item{start}{Bin start position}"
	echo "#'   \item{end}{Binend position}"
	echo "#'   \item{bases}{Percentage of non-N nucleotides acros the full bin}"
	echo "#'   \item{gc}{Percentage of C and G nucleotides \(of non-N nucleotides\)}"
		echo "#'   \item{mappability}{Average mappability of $( echo $KBP | sed -E 's/[PESR]+//' ) with a maximum of 2 mismatches}"
		echo "#'   \item{blacklist}{Percent overlap with High-Signal regions from ExcludeRanges}"
	echo "#'   \item{residual}{Median loess residual from 1000 genomes **PENDING**}"
	echo "#'   \item{use}{Weather bin exceeds minimum thresholds for analysis}"
	echo "#' }"
	echo "#'"
	echo "#' @author RGM"
	echo "#' @source \url{https://github.com/RodrigoGM/QDNAseq.chm13v2}"
	echo "#' @import QDNAseq"
	echo "#' @importFrom Biobase AnnotatedDataFrame"
	echo "\"$1\""
	echo ""
	echo ""
}

for i in $(ls data/ | sed -e 's/.rda//' | sort -k1V) ; do 
	bin_data_documentation $i 
done > R/QDNAseq.chm13v2.R 

```


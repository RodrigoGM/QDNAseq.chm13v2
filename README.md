# QDNAseq.chm13v2: QDNAseq bin annotations for T2T CHM13v2

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R build
status](https://github.com/r-lib/usethis/workflows/R-CMD-check/badge.svg)](https://github.com/r-lib/usethis/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/QDNAseq.chm13v2)](https://CRAN.R-project.org/package=QDNAseq.chm13v2)
<!-- badges: end -->

We provide bin indices to use with QDNAseq at 10, 20, 50, 100, 150, 200, and 500 Kb for the T2T CHM13v2 assembly.  The annotations were created based on the steps from the [QDNAseq vignette](https://bioconductor.org/packages/release/bioc/html/QDNAseq.html), and [QDNAseq.hg38](https://github.com/asntech/QDNAseq.hg38).


**NOTE:  PE100 residuals have now been estimated"

**NOTE:  SR50 bins are pending residual variance estimation with 1000 Genomes Samples. Genomes are currently being downloaded and aligned**

--- 

## Installation

You can install the development version of QDNAseq.chm13v2 like so:

``` r
devtools::install_github("RodrigoGM/QDNAseq.chm13v2", dependencies = TRUE)
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

### 3. 1000 Genomes used for residual estimation
Residual estimation was performed on 38 genomes from the [Simmons Genome Diversity Project](https://www.simonsfoundation.org/simons-genome-diversity-project/) of the [1000 Genomes Project](https://www.internationalgenome.org) via the Data Portal in accord to the [data resuse policy](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/simons_diversity_data/README_Simons_diversity_datareuse_statement.md) (Table 1).  The manifest used to automate the data download is stored in the package, and can be accesed via `data(used.fastq.manifest)`.

For PE100 Genomes were aligned to the `t2t_chm13v2.vd1.fa` assembly using [bwa](https://github.com/lh3/bwa), sorted by [samtools sort](https://github.com/samtools/samtools), and duplicates marked using [Picard MarkDuplicates](https://github.com/broadinstitute/picard), briefly:

```bash
## Main Alignment
bwa mem -aM -t 20 ./t2t_chm13v2.vd1.fa \
	 ${SAMPLE}_1.fastq.gz ${SAMPLE}_2.fastq.gz |
	samtools view -bS - > $sample.chm13v2.bam>


## collate -- sort by read name
samtools collate -@ 20 \
	      -o ${SAMPLE}.chm13v2.cl.bam \
	      ${SAMPLE}.chm13v2.bam

## fixmate
samtools fixmate -m -O BAM -@ 20 \
	      ${SAMPLE}.chm13v2.cl.bam \
	      ${SAMPLE}.chm13v2.fx.bam

## sort bam file
samtools sort -@ 20 \
	      -m 2g \
	      -T ${TMPDIR} \
 	      -o ${SAMPLE}.sorted.bam \
 	      -O BAM \
 	      ${SAMPLE}.chm13v2.fx.bam
		  
## mark duplicate reads
picard MarkDuplicates I=${SAMPLE}.sorted.bam \
	    O=${SAMPLE}.chm13v2.PE.md.bam \
	    M=${SAMPLE}.chm13v2.PE.metrics.txt

## index MD file
samtools index -@ 20 ${SAMPLE}.chm13v2.PE.md.bam

```

#### Table 1.  Samples downloaded for bin residual estimation.
| Sample       | Accession  |
|:-------------|:-----------|
| HG00126      | ERR1025620 |
| HG00128      | ERR1347661 |
| HG00190      | ERR1346534 |
| HG01504      | ERR1025651 |
| HG01600      | ERR1025637 |
| HG01846      | ERR1025638 |
| HG02494      | ERR1395564 |
| HG02724      | ERR1395568 |
| HG02783      | ERR1025663 |
| HG02790      | ERR1025664 |
| HG02943      | ERR1025622 |
| HG03006      | ERR1347669 |
| HG03007      | ERR1347676 |
| HG03085      | ERR1347678 |
| HG03100      | ERR1025621 |
| HGDP00195    | ERR1025649 |
| HGDP00208    | ERR1419159 |
| HGDP00796    | ERR1419128 |
| HGDP00903    | ERR1025606 |
| HGDP01078    | ERR1419130 |
| ALB212       | ERR1395585 |
| AV-21        | ERR1425293 |
| CHI-034      | ERR1347692 |
| Kayseri23827 | ERR1395587 |
| Kor82        | ERR1347707 |
| NA13616      | ERR1347735 |
| NA17377      | ERR1347668 |
| NA18940      | ERR1395570 |
| NA19023      | ERR1347662 |
| NOR111       | ERR1395565 |
| Nesk_22      | ERR1347724 |
| Nesk_25      | ERR1347733 |
| Nlk3         | ERR1347690 |
| TZ-11        | ERR1395617 |
| Utsa21       | ERR1395569 |
| Utsa22       | ERR1395616 |
| ch113        | ERR1025599 |
| mg27         | ERR1025623 |


### 3. Excluded Regions
Exclude regions were obtained from R/Bioconductor [ExcludeRanges](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10126321/) package, [here](https://drive.google.com/drive/folders/1sF9m8Y3eZouTZ3IEEywjs2kfHOWFBSJT).


### 4. Preparing Bins
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

kmers <- c(50, 100) 
names(kmers) <- paste0("PE", kmers)
names(kmers)[1] <- "SR50"

t2t.mappability <- "../Mappability/"

( bam.files= list(
    "SR50" = NULL,
    "PE100" = list.files(path = "data-raw/bwa_out/PE100/", pattern = "PE.md.bam$",
                         full.names = TRUE, recursive = TRUE))
)


## loop through PE100 
for(k in names(kmers[2])) {
    ## loop through bin.sizes
    for(binsize in bin.sizes) {

        message(paste(date(), ": starting", binsize, "kbp run"))
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

        message(t2t.mp.bigwig)
        ## There was a bug in the calculateMappability, thus
        ## computing directly 
		
		## comand line construction
        cmd <- paste(bigWigAvgOB, t2t.mp.bigwig, bed.file,
                     "stdout | cut -f 1,5 ")
					 
		## run cmd, internalize ouptut, and convert to data.frame
        mappability <- system(cmd, intern = TRUE)
        mappability <- as.data.frame(
            do.call(rbind,
                    strsplit(mappability, "\t"))) %>%
            transform(row.names = V1, V1 = NULL) %>%
            rename("mappability" = "V2")
			
        ## coerce output to numberic, and convert to percent
        bins$mappability <- as.numeric(
            mappability[rownames(bins), "mappability"]) * 100

	    ## coerce GC content to numeric
        bins$gc <- as.numeric(bins$gc)

        ## estimate % overap to exclude list
        bins$blacklist <- calculateBlacklist(
            bins,
            bedFiles = t2t.exclude.file)

        ## make residual column
        bins$residual <- NA

        ## make use column
        bins$use <- bins$chromosome %in% as.character(1:22) & bins$bases > 0

        ## Count bins across 1000 Genomes bam files
        tg <- binReadCounts(bins,
                            bamfiles = bam.files[[k]],
                            cache=TRUE,
                            isPaired = TRUE,
                            pairedEnds = TRUE)
							
	    ## estimate residuals
        bins$residual <- iterateResiduals(tg)
		
		## convert to annotated data frame
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
                    "Median loess residual from 1000 Genomes (PE 100mers alignment)",
                    "Whether the bin should be used in subsequent analysis steps"),
                row.names = colnames(bins)))

 	    ## metadata
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

        message(paste(date(), ": end of", binsize, "kbp run"))

    }
}

```

## -- tidy up -- 
The following bash routine was used to change object names from `bins` to match the data, and to document each data set

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


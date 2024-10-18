## code to generate bins on T2T CHM13v2
## use R/R-4.1.2
## module load R/R-4.1.2
if(R.home() != "/opt/common/CentOS_7/R/R-4.1.2/lib64/R") {
    stop("run `module load R/R-4.1.2`")
}


## libraries
library(Biobase)
library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)
library(QDNAseq)
library(future)
library(tidyverse)

#set virtual mem
options(future.globals.maxSize= 120*1023^3) #8912896000)

# change the current plan to access parallelization
future::plan("multicore", workers = 11)

## ---- prepare bins at various sizes

bin.sizes <- c(1000, 500, 200, 150, 100, 50, 30, 20, 15, 10, 5, 1)
names(bin.sizes) = paste0(bin.sizes, "kbp")

kmers <- c(50, 100) ## 1000 genomes does not have PE 150 data
names(kmers) <- paste0("PE", kmers)
names(kmers)[1] <- "SR50"

t2t.mappability <- "../Mappability/"

( bam.files= list(
    "SR50" = NULL,
    "PE100" = list.files(path = "data-raw/bwa_out/PE100/", pattern = "PE.md.bam$",
                         full.names = TRUE, recursive = TRUE))
)

## loop through SR50, PE100 and PE150 
for(k in names(kmers[2])) {
    ## loop through bin.sizes
    for(binsize in bin.sizes) {
        ##
        message(paste(date(), ": starting", binsize, "kbp run"))
        ## create bins
        bins <- createBins(bsgenome = BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0,
                           binSize = binsize, ignoreMitochondria = FALSE)
        ##        
        ## write out bed file
        bed.file <- paste0("chm13v2.", binsize, "kbp.bed")
        write.table(cbind(bins[,1:3], rownames(bins)),
                    file = bed.file,
                    sep = "\t",
                    row.names = FALSE, col.names = FALSE, quote = FALSE)
        ##
        ## mappability files
        t2t.mp.bigwig <- file.path(
            t2t.mappability,
            paste0("t2t_genmap_k", kmers[k],
                   "_E2/chm13v2.vd1.k", kmers[k], "_E2.genmap.bw"))
        ##
        ## exclude list file
        t2t.exclude.file <- file.path(t2t.mappability,
                                      "T2T.excluderanges.bed")
        bigWigAvgOB <- file.path("bigWigAverageOverBed")
        ##
        message(t2t.mp.bigwig)
        ## there was a bug in the calculateMappability, thus
        ## computing directly and reading in file
        ## compute mapability per bin
        ## mappability.file <- paste0("chm13v2.", binsize, "kbp.", k, ".bed")
        cmd <- paste(bigWigAvgOB, t2t.mp.bigwig, bed.file,
                     "stdout | cut -f 1,5 ")
        mappability <- system(cmd, intern = TRUE)
        mappability <- as.data.frame(
            do.call(rbind,
                    strsplit(mappability, "\t"))) %>%
            transform(row.names = V1, V1 = NULL) %>%
            rename("mappability" = "V2")                    
        bins$mappability <- as.numeric(
            mappability[rownames(bins), "mappability"]) * 100
        ##
        bins$gc <- as.numeric(bins$gc)
        ##
        ## estimate % overap to exclude list
        bins$blacklist <- calculateBlacklist(
            bins,
            bedFiles = t2t.exclude.file)
        ##
        ## make residual column
        bins$residual <- NA
        ##
        ## make use column
        bins$use <- bins$chromosome %in% as.character(1:22) & bins$bases > 0
        ##
        ## Count bins across 1000 Genomes bam files
        tg <- binReadCounts(bins,
                            bamfiles = bam.files[[k]],
                            cache=TRUE,
                            isPaired = TRUE,
                            pairedEnds = TRUE)
        ##
        bins$residual <- iterateResiduals(tg)
        ##
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
        ##
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
        ##
        attr(bins, "QDNAseq") <- QDNAseqInfo
        save(bins, file=file.path("data",
                                  paste0("chm13v2.", binsize, "kbp.", k,".rda")),
             compress='xz')
        ##
        message(paste(date(), ": end of", binsize, "kbp run"))
                ##
    }
}


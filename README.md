# **quaffle**: Quantifying Usage of Alternative First and Last Exons

## How to install
```r
# install.packages("devtools")
devtools::install_github("fursham-h/quafle")
```

## What you need
1. Reference annotation (GTF format)
2. Directory containing BAM files (and its indices)

## Usage
```r
# load quafle into environment
library(quafle)

# Step 1: build a database of alternative first and last exons
gtf <- "/path/to/gtf"
afl.db <- buildAFL(gtf)

# Step 2: Quantify PSI values for each sample
bams <- "/path/to/bamdir"  # This directory should contain all bam files and its indices
afl.psi <- quantAFL(afl.db, bams)
# a tsv file with PSI values will be saved in working directory

# Step 3: Compare exon usage between two samples
## create a character vector with names of replicates
## Note: these names have to match the name of BAM files, without .bam extension
a <- c("sample_A_1", "sample_A_2")
b <- c("sample_B_1", "sample_B_2")
afl.diff <- diffAFL(afl.psi, a, b)

# Optional step: Retrieving exon sequence
## Load genome
library(BSgenome.Mmusculus.UCSC.mm10)

# run function
seq <- getAFLseq(afl.diff, Mmusculus)



```



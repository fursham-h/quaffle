# **Quaffle**: Quantifying Usage of Alternative First and Last Exons

## How to install
```r
# install.packages("devtools")
devtools::install_github("fursham-h/quaffle")
```

## What you need
1. Reference annotation (GTF format)
2. Directory containing BAM files (and its indices)
3. DataFrame containing samples and its grouping

## Usage
```r
# load Quaffle into environment
library(quaffle)

# Step 0: Prepare sample DataFrame and variables

## Sample names should not contain .bam extension
## Set names of rows to reflect sample names
sampledata <- data.frame(samples = c("SampleA1","SampleA2","SampleB1", "SampleB2"),
                      group = c("A","A","B","B")    
                      )
rownames(sampledata) <- sampledata$samples

# set GTF directory and BAM directory
gtf <- "/path/to/gtf"
bams <- "/path/to/bamdir"  # This directory should contain all bam files and its indices

# Step 1: Create quaffle object
## This creates ranges of all AFL events,
## count reads mapping to AFL exons,
## and quantify PSI values
qobj <- createQuaffle(gtf, 
                      colData = sampledata, 
                      bamdir = bams)


# Step 2: Compare exon usage between two samples
## Provide a character vector with the following information:
### 1. Name of factor for sample A
### 2. Name of factor for sample B
### 3. Variable name in colData containing sample groupings

contrast.vector <- c("A","B","group")

# run differential analysis
afl.diff <- diff(qobj, contrast = contrast.vector)

# Optional step: Retrieving exon sequence
## Load genome
### The easiest way is to use genome from BSgenome package
library(BSgenome.Mmusculus.UCSC.mm10)

# run function
seq <- getAFLseq(afl.diff, Mmusculus)

# sequences can be exported as fasta files by providing
# pathname to `dir` argument
seq <- getAFLseq(afl.diff, Mmusculus, dir = "output")


```



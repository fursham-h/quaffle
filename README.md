# **Quaffle**: Quantifying Usage of Alternative First and Last Exons

## How to install
```r
# install.packages("devtools")
devtools::install_github("fursham-h/quaffle")
```

## What you need
1. Reference annotation (GTF format)
2. Directory containing BAM files (and its indices)

## Usage
Below is an example workflow, using sample datasets that come with Quaffle
```r
# load Quaffle into environment
library(quaffle)

# set GTF directory and BAM directory
gtf <- system.file("extdata/wtap.gtf", package = "quaffle")
bams <- system.file("extdata/bams", package = "quaffle")

# Step 1: Run Quaffle
## This creates ranges of all AFL events,
## count reads mapping to AFL exons,
## and quantify PSI values
se <- Quaffle(bams, gtf)


# Step 2: Run differential usage (Bludger)
groupings <- rep(c("A","B"), each = 3)
results <- Bludger(se, groupings)

# Optional step: Retrieving AFL exon sequence
## Load genome
### The easiest way is to use genome from BSgenome package
library(BSgenome.Mmusculus.UCSC.mm10)

# run function
seq <- Snitch(se, Mmusculus)

```



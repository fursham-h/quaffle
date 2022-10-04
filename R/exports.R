#' Run Quaffle
#'
#' @description `Quaffle` will identify all alternative first and last (AFL) exons
#' and quantify the length-normalized read counts for each sample. This function
#' require a reference GTF file from which an AFL database is built upon, and
#' a directory containing aligned read files (BAMs). `Quaffle` outputs a
#' RangedSummarizedExperiment object that stores AFL coordinates, normalized
#' read counts and sample metadata
#'
#' @param bamdir Directory containing BAM files and indices
#' @param gtf Path to reference GTF file
#' @param colData Optional: A dataframe containing names of samples as rownames
#' and other sampel information
#'
#' @return RangedSummarizedExperiment object, with 2 assays; inc and total normalized
#' read counts. Number of columns in object equals to number of samples/BAM files
#' and number of rows equals to number of identified alternative First and Last
#' exons.
#' @export
#'
#' @examples
#' library(quaffle)
#' gtf <- system.file("extdata/wtap.gtf", package = "quaffle")
#' bams <- system.file("extdata/bams", package = "quaffle")
#'
#' se <- Quaffle(bams, gtf)
Quaffle <- function(bamdir, gtf, colData=NULL,
                    pairedEnd = FALSE, nthreads = 1, strandness = 0){

    # run input checks
    mandargs <- c("gtf", "bamdir")
    passed <- names(as.list(match.call())[-1])
    if (any(!mandargs %in% passed)) {
        rlang::abort(paste(
            "missing values for",
            paste(setdiff(mandargs, passed), collapse = ", ")
        ))
    }

    if(!file.exists(gtf)){
        rlang::abort("GTF annotation does not exist")
    }
    if(!dir.exists(bamdir)){
        rlang::abort("BAM directory does not exist")
    }
    if(!is.null(colData)){
        bams <- paste0(rownames(colData),".bam")
        if(!all(bams %in% list.files(bamdir))){
            rlang::abort(sprintf("BAM files not found at %s directory", bamdir))
        }
    } else {
        bams <- stringr::str_remove(list.files(bamdir, pattern = ".bam$"), ".bam")
        if(length(bams) == 0){
            rlang::abort(sprintf("No BAM files found at %s directory", bamdir))
        }
        colData <- data.frame(samples = bams,
                              row.names = bams)
    }

    # create se object
    ranges <- .buildAFL(gtf)
    counts <- .countAFL(ranges, colData, bamdir, pairedEnd, nthreads, strandness)
    x <- SummarizedExperiment::SummarizedExperiment(
        assays = counts,
        rowRanges = ranges,
        colData = colData
    )

    return(x)
}


#' Run Double Binomial GLM differential analyses
#'
#' @description `Bludger` compares the usage of each AFL event between two groups,
#' using a double binomial GLM model.
#'
#' @param se RangedSummarizedExperiment object from `Quaffle` output
#' @param group Can be name of variable from `colData(se)` that contain groupings
#' of samples. Can also be a vector (length equal to number of samples) that
#' describes the sample groupings
#' @param contrast numeric vector of length 2 specifying which levels of
#' the "groups" factor should be compared.
#'
#' @return a dataframe describing the output of differential AFL usage with
#' the folloinwing column names:
#'  \itemize{
#'  \item{gene_id}{ID of parent gene}
#'  \item{gene_name}{Name of parent gene}
#'  \item{effectivecoord}{Coordinate of AFL event}
#'  \item{type}{Type of event; Alternative First (AF) or Last (AL)}
#'  \item{MLE_*}{Maximum likelihood estimate of event usage in each sample}
#'  \item{deltaMLE}{Change in event usage between comparisons}
#'  \item{pval}{Unadjusted p-value of differential analysis}
#'  \item{adj_pval}{Adjusted p-value of differential analysis}
#'  \item{MeanTotCt_*}{Mean total read count in each sample}
#' }
#' @export
#'
#' @examples
#' #' library(quaffle)
#' gtf <- system.file("extdata/wtap.gtf", package = "quaffle")
#' bams <- system.file("extdata/bams", package = "quaffle")
#'
#' se <- Quaffle(bams, gtf)
#'
#' Bludger(se, rep(c("A","B"), each = 3))
#'
Bludger <- function(se, group, contrast = c(1,2)){

    # run input checks
    mandargs <- c("se", "group")
    passed <- names(as.list(match.call())[-1])
    if (any(!mandargs %in% passed)) {
        rlang::abort(paste(
            "missing values for",
            paste(setdiff(mandargs, passed), collapse = ", ")
        ))
    }

    # check object type
    if(!is(se, "RangedSummarizedExperiment") | any(names(se@assays) != c("inc", "total"))){
        rlang::abort(str_glue("Unrecognized input object"))
    }

    #check groupings
    if(typeof(group)=="character" & length(group) == 1 ){
        if(group %in% colnames(se@colData)){
            group <- se[[group]]
        } else {
            rlang::abort(str_glue("`group` variable not found in SummarizedExperiment object"))
        }
    } else if(length(group) != ncol(se)){
        rlang::abort(str_glue("`group` length ({length(group)}) not equal to number of samples ({ncol(se)})"))
    }


    dge <- DoubleExpSeq::DBGLM1(se@assays@data$inc, se@assays@data$total, groups = group,
                                contrast = contrast)
    meta_to_join <- se@rowRanges %>%
        as.data.frame() %>%
        dplyr::select(effectivecoord, gene_id, gene_name, type)
    out <- dge$All %>%
        as.data.frame() %>%
        cbind(meta_to_join) %>%
        dplyr::mutate(deltaPSI = .[[1]] - .[[2]]) %>%
        dplyr::select(gene_id, gene_name, effectivecoord, type,
                      starts_with("MLE"), deltaPSI, pval = pVal,
                      adj_pval = Adj.pVal, starts_with("Mean"))
    rownames(out) <- NULL
    out

}


#' Get sequence
#'
#' @param se RangedSummarizedExperiment object from `Quaffle` output
#' @param fasta BSGenome object or DNAStringList object containing genome sequence
#' @param subset A character vector containing IDs to subset. IDs can be name of
#' genes, ID of genes or AFL coordinates
#' @param type Type of sequence to return. Can be the start of exon ("start"),
#' end of exon ("end") or entire exon ("exon")
#' @param exontype Type of event to return, Can be "all", "AF" or "AL"
#' @param upstream Length of upstream padding sequence
#' @param downstream Length of downstream padding sequence
#' @param outdir Path to output directory (Default: only return DNAStringSet object)
#' @param filename Name of FASTA file to create
#'
#' @return DNAStringSet object with DNA sequences of selected events
#' @export
Snitch <- function(se, fasta, subset = NULL,
                   type = c("start","end","exon"),
                   exontype = c("all", "AF", "AL"),
                   upstream = 50, downstream = 50,
                   outdir = NULL, filename = "Snitch.fasta"){

    # Checks
    # catch missing args
    mandargs <- c("se", "fasta")
    passed <- names(as.list(match.call())[-1])
    if (any(!mandargs %in% passed)) {
        rlang::abort(paste(
            "missing values for",
            paste(setdiff(mandargs, passed), collapse = ", ")
        ))
    }
    if(any(!type %in% c("start", "end", "exon"))){
        rlang::abort(sprintf("Unknown type `%s`", type))
    }
    if(!is.null(outdir)){
        if(!dir.exists(outdir)){
            dir.create(outdir)
        }
    }

    type <- type[1]
    exontype <- exontype[1]

    # get GRanges object
    coord <- GenomicRanges::GRanges(se@rowRanges$effectivecoord,
                                    strand = strand(se@rowRanges))
    coord$effectivecoord <- se@rowRanges$effectivecoord
    coord$name <- paste0(se@rowRanges$gene_name, "_",
                         se@rowRanges$effectivecoord,"_",
                         se@rowRanges$type, "_")

    if(toupper(exontype) != "ALL"){
        coord <- coord[se@rowRanges$type==exontype]
    }

    if(!is.null(subset)){
        coord <- coord[coord$effectivecoord %in% .getcoords(se@rowRanges, subset)]
    }


    # get desired coord
    if(type == "start"){
        coord <- GenomicRanges::resize(coord, 1)
    } else if(type == "end"){
        coord <- GenomicRanges::resize(coord, 1, fix = "end")
    }
    coord <- GenomicRanges::resize(coord,
                                   BiocGenerics::width(coord)+upstream,
                                   fix = "end")
    coord <- GenomicRanges::resize(coord,
                                   BiocGenerics::width(coord)+downstream,
                                   fix = "start")

    seq <- Biostrings::getSeq(fasta, coord)
    names(seq) <- coord$name

    if(!is.null(outdir)){
        Biostrings::writeXStringSet(x = seq,
                                    filepath = file.path(outdir, filename))
    }
    return(seq)
}































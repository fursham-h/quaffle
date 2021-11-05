#' Create Quaffle Object
#'
#' @param gtf Path to GTF annotation
#' @param colData DataFrame containing sample information
#' @param bamdir Path to directory containing BAM files. It is preferable to have
#' bam indices (.bam.bai) in the same directory.
#' @param min_read Minimum number of aligned read on AFL exons
#'
#' @return Quaffle object
#' @export
#'
#'
createQuaffle <- function(gtf,
                          colData,
                          bamdir,
                          min_read=5){


    mandargs <- c("gtf", "colData", "bamdir")
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
    bams <- paste0(rownames(colData),".bam")
    if(!all(bams %in% list.files(bamdir))){
        rlang::abort(sprintf("BAM files not found at %s directory", bamdir))
    }

    # create object
    x <- methods::new("QuaffleObject")
    colData(x) <- colData
    x <- build(x, gtf)
    x <- count(x, bamdir)
    x <- quant(x, min_read)
    return(x)
}

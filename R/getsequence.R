#' Retrieve feature sequences
#'
#' @param diff
#' DataFrame containing differential splicing analysis, output from diff() function
#' @param fasta
#' XString object containing genome sequence
#' @param dir
#' Output directory for FASTA files (Optional)
#' @param type
#' Type of sequence to output. Can be "start",
#' "end" or "exon" which will output the 5' border,
#' 3' end or the entire exon sequence respectively.
#' @param exontype
#' Type of exon to output. Can be "all" (default),
#' "AF" or "AL" exons.
#' @param upstream
#' Amount of padding to extract for upstream sequence (Default: 50bp)
#' @param downstream
#' Amount of padding to extract for downstream sequence (Default: 50bp)
#' @param p.value
#' Threshold to select for statistically significant events. By default,
#' it will use a value of 0.05 on adjusted p-values (adj_p_val). See
#' `use.unadjusted` parameter to select on unadjusted p-values.
#' @param deltapsi
#' Threshold to select for statistically significant events. By default,
#' will select events with 10% change in either samples (0.1)
#' @param use.undajusted
#' Boolean value on whether to use unadjusted p-values
#'
#' @return
#' @export
#'
#' @examples
getAFLseq <- function(diff, fasta,
                      dir = NULL,
                      type = c("start","end","exon"),
                      exontype = c("all", "AF", "AL"),
                      upstream = 50,
                      downstream = 50,
                      p.value = 0.05,
                      deltapsi = 0.1,
                      use.undajusted = F){


    p_val <- adj_p_val <- deltaPSI <- NULL

    # Checks
    # catch missing args
    mandargs <- c("diff", "fasta")
    passed <- names(as.list(match.call())[-1])
    if (any(!mandargs %in% passed)) {
        rlang::abort(paste(
            "missing values for",
            paste(setdiff(mandargs, passed), collapse = ", ")
        ))
    }
    if(!type %in% c("start", "end", "exon")){
        rlang::abort(sprintf("Unknown type `%s`", type))
    }
    if(!is.null(dir)){
        if(!dir.exists(dir)){
            dir.create(dir)
        }
    }

    type <- type[1]
    exontype <- exontype[1]

    # prepare diff
    if(use.undajusted){
        diff <- dplyr::filter(diff, p_val < p.value)
    } else {
        diff <- dplyr::filter(diff, adj_p_val < p.value)
    }
    if(exontype!="all"){
        diff <- dplyr::filter(diff, type == exontype)
    }
    diff <- diff %>%
        dplyr::filter(abs(deltaPSI) > deltapsi) %>%
        dplyr::mutate(sample = ifelse(deltaPSI > 0, "B", "A"))


    coord <- GenomicRanges::GRanges(diff$effectivecoord,
                                    strand = diff$strand)

    # get desired coord
    if(type == "start"){
        coord <- GenomicRanges::resize(coord, 1)
        c
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
    names(seq) <- paste0(diff$gene_name, "_",diff$effectivecoord)

    A <- seq[which(diff$sample == "A")]
    B <- seq[which(diff$sample == "B")]



    if(!is.null(dir)){
        Biostrings::writeXStringSet(x = A,
                                    filepath = paste0(dir, "/A.fasta"))
        Biostrings::writeXStringSet(x = B,
                                    filepath = paste0(dir, "/B.fasta"))
    }
    return(list(A = A, B = B))

}

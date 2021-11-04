getAFLseq <- function(diff, fasta,
                      dir = NULL,
                      type = "start",
                      exontype = "all",
                      upstream = 50,
                      downstream = 50,
                      p.value = 0.05,
                      deltapsi = 0.1,
                      use.undajusted = F){

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
    } else {
        return(list(A = A, B = B))
    }

}

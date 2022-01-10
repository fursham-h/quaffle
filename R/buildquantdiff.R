#' @title Build list of alternate first and last exons
#'
#' @param x
#' Path to GTF annotation
#'
#' @return
#' GenomicRanges object with coordinates of alternative first and last exons
#' extracted from the input annotation
#'
#'
#' @importFrom dplyr %>%
#' @importFrom S4Vectors DataFrame
#'
.buildAFL <- function(x) {

    end <- transcript_id <- strand <- start <- afl.seqnames <- afl.start <- NULL
    afl.end <- afl.strand <- gene_id <- gene_name <- type <- coding <- NULL

    # import transcriptome
    rlang::inform("Reading GTF input")
    x <- rtracklayer::import(x)
    # check for CDS entries
    assertthat::assert_that(
        "CDS" %in% unique(x$type),
        msg = "Input GTF is missing CDS information")
    x.cds <- x[x$type == "CDS"]
    x <- x[x$type == "exon"]

    # get list of first exons from each transcript
    rlang::inform("Building AFL database")
    all.start <- S4Vectors::split(x, ~transcript_id) %>%
        range() %>%
        GenomicRanges::resize(1) %>%
        unlist() %>%
        GenomicRanges::reduce()

    cds.afl <- x.cds %>% as.data.frame() %>%
        dplyr::group_by(transcript_id) %>%
        dplyr::arrange(ifelse(strand == "-", dplyr::desc(start), start)) %>%
        dplyr::filter(dplyr::row_number() == 1 | dplyr::row_number() == dplyr::n()) %>%
        dplyr::ungroup() %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T) %>%
        GenomicRanges::reduce()

    # get a list of all first and last exons from annotation
    afl <- x %>% as.data.frame() %>%
        dplyr::group_by(transcript_id) %>%
        dplyr::arrange(ifelse(strand == "-", dplyr::desc(start), start)) %>%
        dplyr::filter(dplyr::row_number() == 1 | dplyr::row_number() == dplyr::n()) %>%
        dplyr::ungroup() %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T) %>%
        GenomicRanges::reduce()

    # annotate AF and AL
    afl$type <- ifelse(IRanges::overlapsAny(afl, all.start, type = "start") |
                          IRanges::overlapsAny(afl, all.start, type = "end"),
                      "AF", "AL")
    # annotate CDS
    afl$coding <- ifelse(IRanges::overlapsAny(afl, cds.afl, type = "start") |
                          IRanges::overlapsAny(afl, cds.afl, type = "end"),
                      "TRUE", "FALSE")

    # annotate gene_id and gene_name of each exons
    afl.labelled <- IRanges::mergeByOverlaps(afl, x) %>%
        as.data.frame() %>%
        dplyr::select(seqnames = afl.seqnames, start = afl.start, end = afl.end,
                      strand = afl.strand, gene_id, gene_name, type, coding) %>%
        dplyr::distinct() %>%
        dplyr::group_by(gene_id, type) %>%
        dplyr::filter(dplyr::n() > 1) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(ifelse(strand == "-", dplyr::desc(start), start)) %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
    afl.labelled$effectivecoord <- afl.labelled$coord <- paste0(GenomeInfoDb::seqnames(afl.labelled),
                                          ":",
                                          BiocGenerics::start(afl.labelled),
                                          "-",
                                          BiocGenerics::end(afl.labelled))

    rlang::inform("Extracting effective coordinates")
    int.exons <- x %>% as.data.frame() %>%
        dplyr::group_by(transcript_id) %>%
        dplyr::arrange(ifelse(strand == "-", dplyr::desc(start), start)) %>%
        dplyr::filter(dplyr::row_number() != 1 & dplyr::row_number() != dplyr::n()) %>%
        dplyr::ungroup() %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T) %>%
        GenomicRanges::reduce()

    internal <- IRanges::findOverlapPairs(afl.labelled, int.exons)
    exclude <- BiocGenerics::end(internal@second) < BiocGenerics::end(internal@first) & BiocGenerics::start(internal@second) > BiocGenerics::start(internal@first)
    internal <- internal[!exclude]
    changed <- IRanges::findOverlaps(internal@first,
                            afl.labelled, type = "equal",
                            select = "first")
    internal.setdiff <- GenomicRanges::psetdiff(internal)
    afl.labelled$effectivecoord[changed] <- paste0(GenomeInfoDb::seqnames(internal.setdiff),
                                          ":",
                                          BiocGenerics::start(internal.setdiff),
                                          "-",
                                          BiocGenerics::end(internal.setdiff))

    return(afl.labelled)
}


#' Quantify reads mapping to alternate first and last exons
#'
#' @param x Quaffle object
#' @param dir Path to directory containing BAM files. It is preferable to have
#' bam indices (.bam.bai) in the same directory.
#'
#' @return
#' Data-frame containing read count and PSI for each AFL
#'
#'
.countAFL <- function(x, dir) {

    count <- width <- gene_id <- type <- norm_count <- totalnormcount <- NULL
    strand <-  start <-  seqnames <-  coding <-  PSI <- NULL

    # Checks
    # catch missing args

    # retrieve list of bam files
    db <- x@rowRanges
    dir <- stringr::str_remove(dir, "/$")
    bams <- paste0(rownames(x@colData),".bam")
    if(!all(bams %in% list.files(dir))){
        rlang::abort(sprintf("BAM files not found at %s directory", dir))
    }

    # Prepare parameters for scanBam
    which <- db
    what <- c("rname", "strand", "pos")
    param <- Rsamtools::ScanBamParam(which=which, what=what)

    # loop to:
    ## 1) read BAM files and convert it to a GR object
    ## 2) calculate number of reads mapping to each AFL
    rlang::inform("Quantifying AFL reads")
    out <- BiocParallel::bplapply(bams, function(x){
        xpath <- paste0(dir, "/", x)


        bam <- Rsamtools::scanBam(xpath, param = param)
        bam <- unname(bam)
        elts <- stats::setNames(Rsamtools::bamWhat(param),
                                Rsamtools::bamWhat(param))
        lst <- lapply(elts, function(elt) .bamunlist(lapply(bam, "[[", elt)))
        bam.gr <- do.call("DataFrame", lst) %>%
            GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "rname",
                                                    start.field = "pos", end.field = "pos")

        rlang::inform(sprintf("\n\tSuccessfully read %s file", x))
        return(IRanges::countOverlaps(db, bam.gr))
    })
    out <- suppressMessages(out %>%
        dplyr::bind_cols()) %>%
        as.matrix()

    # Correct sample names
    colnames(out) <- stringr::str_remove(bams, ".bam$")
    rownames(out) <- db$coord
    out



}

#' Quantify PSI metrics
#'
#' @param db rowRanges slot from Quaffle object
#' @param out counts slot from Quaffle object
#'
#' @return Matrix containing PSI values
.quantAFL <- function(db, out){

    width <- gene_id <- type <- norm_count <- totalnormcount <- NULL
    strand <- start <- coord <- PSI <- NULL

    # clean data and calculate normalized counts and PSI
    out.comb <- as.data.frame(db) %>%
        dplyr::bind_cols(as.data.frame(out)) %>%
        tidyr::gather("sample", "count", colnames(out)) %>%
        #dplyr::mutate(state = ifelse(count < min_read, "LOW", "OK")) %>%
        dplyr::mutate(norm_count = count/width) %>%
        dplyr::group_by(gene_id, sample, type) %>%
        dplyr::mutate(totalnormcount = sum(norm_count)) %>%
        dplyr::mutate(PSI = signif(norm_count/totalnormcount, digits = 3)) %>%
        tidyr::replace_na(list(PSI = 0)) %>%
        dplyr::select(-norm_count, -totalnormcount) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(gene_id) %>%
        dplyr::arrange(ifelse(strand == "-", dplyr::desc(start), start),
                       type, sample) %>%
        dplyr::ungroup()


    psi <- out.comb %>%
        dplyr::select(coord, sample, PSI) %>%
        dplyr::distinct(coord, sample, .keep_all = T) %>%
        tidyr::spread(sample, PSI) %>%
        tibble::column_to_rownames("coord") %>%
        as.matrix()


    return(psi)
}

.diff <- function(x, colData, state, comparison, min_samples, adjust.methods){

    id <- psilist <- B.PSI <- A.PSI <- p_val <- deltaPSI <- adj_p_val <- NULL
    grouping <- comparison[[3]]
    A <- comparison[[1]]
    B <- comparison[[2]]
    state <- state[rownames(x),]
    state <- state %>%
        as.data.frame() %>%
        #dplyr::mutate(id = rownames(x)) %>%
        tibble::rownames_to_column("id") %>%
        tidyr::gather("samples","state", -id)

    x <- .posterior(x) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("id") %>%
        tidyr::gather("samples","psi", -id)



    colData %>%
        dplyr::left_join(x, by = "samples") %>%
        dplyr::left_join(state, by = c("samples", "id")) %>%
        dplyr::group_by(id, !!!rlang::syms(grouping)) %>%
        dplyr::filter(sum(state == "OK") >= min_samples) %>%
        dplyr::summarise(psilist = list(psi)) %>%
        tidyr::spread(grouping, psilist) %>%
        dplyr::filter(!S4Vectors::isEmpty(!!!rlang::syms(A)) & !S4Vectors::isEmpty(!!!rlang::syms(B))) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(p_val = this.t.test(!!!rlang::syms(A),!!!rlang::syms(B))) %>%
        dplyr::mutate(A.PSI = mean(!!!rlang::syms(A)),
               B.PSI = mean(!!!rlang::syms(B)),
               deltaPSI = B.PSI-A.PSI) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(adj_p_val = stats::p.adjust(p_val, method = adjust.methods)) %>%
        dplyr::mutate(A.PSI = signif(A.PSI, digits = 3)) %>%
        dplyr::mutate(B.PSI = signif(B.PSI, digits = 3)) %>%
        dplyr::mutate(deltaPSI = signif(deltaPSI, digits = 3)) %>%
        dplyr::mutate(p_val = signif(p_val, digits = 3)) %>%
        dplyr::mutate(adj_p_val = signif(adj_p_val, digits = 3)) %>%
        dplyr::select(coord = id, A.PSI, B.PSI, deltaPSI, p_val, adj_p_val)


}

this.t.test <- function(...){
    obj <- tryCatch({
             Bolstad::bayes.t.test(...)$p.value},
             error=function(cond) return(1))
    return(obj)
}


.bamunlist <- function (x) {
    x1 <- x[[1L]]
    if (is.factor(x1)) {
        structure(unlist(x), class = "factor", levels = levels(x1))
    } else {
            do.call(c, x)
    }
}

#' Build list of alternate first and last exons
#'
#' @description
#'
#' @param x
#' Path to GTF annotation
#'
#' @return
#' GenomicRanges object with coordinates of alternative first and last exons
#' extracted from the input annotation
#'
#' @export
#'
#' @examples
#' gtf <- system.file("extdata", "wtap.gtf", package = "AFLanalyze")
#' buildAFL(gtf)
#'
buildAFL <- function(x) {
    # import transcriptome
    rlang::inform("Reading GTF input")
    x <- rtracklayer::import(x)
    x <- x[x$type == "exon"]

    # get list of first exons from each transcript
    rlang::inform("Building AFL database")
    all.start <- S4Vectors::split(x, ~transcript_id) %>%
        range() %>%
        GenomicRanges::resize(1) %>%
        unlist() %>%
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
    afl$pos <- ifelse(IRanges::overlapsAny(afl, all.start, type = "start") |
                          IRanges::overlapsAny(afl, all.start, type = "end"),
                      "AF", "AL")

    # annotate gene_id and gene_name of each exons
    afl.labelled <- IRanges::mergeByOverlaps(afl, x) %>%
        as.data.frame() %>%
        dplyr::select(seqnames = afl.seqnames, start = afl.start, end = afl.end,
                      strand = afl.strand, pos, gene_id, gene_name) %>%
        dplyr::distinct() %>%
        dplyr::group_by(gene_id, pos) %>%
        dplyr::filter(dplyr::n() > 1) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(ifelse(strand == "-", dplyr::desc(start), start)) %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)

    return(afl.labelled)
}


#' Quantify reads mapping to alternate first and last exons
#'
#' @param db Database of AFL from buildAFL output
#' @param dir Path to directory containing BAM files. It is preferable to have
#' bam indices (.bam.bai) in the same directory.
#'
#' @return
#' Data-frame containing read count and PSI for each AFL
#'
#' @export
#'
#' @examples
quantAFL <- function(db, dir, min_read = 5) {

    # Checks
    # catch missing args
    mandargs <- c("db", "dir")
    passed <- names(as.list(match.call())[-1])
    if (any(!mandargs %in% passed)) {
        rlang::abort(paste(
            "missing values for",
            paste(setdiff(mandargs, passed), collapse = ", ")
        ))
    }
    argnames <- as.character(match.call())[-1]
    assertthat::assert_that(assertthat::is.dir(dir))
    assertthat::assert_that(is(db, "GRanges"), msg = sprintf("`%s` is not a GRanges object", argnames[1]))

    # retrieve list of bam files
    dir <- stringr::str_remove(dir, "/$")
    bams <- list.files(dir, pattern = ".bam$", full.names = F)

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

        IRanges::countOverlaps(db, bam.gr)
        rlang::inform(sprintf("\n\tSuccessfully read %s file", x))
    })
    out <- suppressMessages(out %>%
        dplyr::bind_cols())

    # Correct sample names
    names(out) <- stringr::str_remove(bams, ".bam$")

    # clean data and calculate normalized counts and PSI
    out.comb <- as.data.frame(db) %>%
        dplyr::bind_cols(out) %>%
        tidyr::gather("sample", "count", names(out)) %>%
        dplyr::mutate(count = ifelse(count < min_read, 0, count)) %>%
        dplyr::mutate(norm_count = count/width) %>%
        dplyr::group_by(gene_id, sample, pos) %>%
        dplyr::mutate(totalnormcount = sum(norm_count)) %>%
        dplyr::mutate(PSI = signif(norm_count/totalnormcount, digits = 3)) %>%
        tidyr::replace_na(list(PSI = 0)) %>%
        dplyr::select(-norm_count, -totalnormcount)

    out.comb %>%
        dplyr::select(seqnames:gene_name, sample,PSI) %>%
        tidyr::spread(sample, PSI) %>%
        write.table("QuaflePSI.tsv", sep = '\t', quote = F, row.names = F)
    return(out.comb)

}


#' Differential AFL analysis
#'
#' @param psi Dataframe with AFL PSI values from quantAFL output
#' @param a Character list containing sample names from group "a". Sample names
#' have to be identical to the name of the BAM files (without .BAM extension)
#' @param b Character list containing sample names from group "b". Sample names
#' have to be identical to the name of the BAM files (without .BAM extension)
#'
#' @return
#' @export
#'
#' @examples
diffAFL <- function(psi, a, b) {

    # Checks
    # catch missing args
    mandargs <- c("psi", "a", "b")
    passed <- names(as.list(match.call())[-1])
    if (any(!mandargs %in% passed)) {
        rlang::abort(paste(
            "missing values for",
            paste(setdiff(mandargs, passed), collapse = ", ")
        ))
    }
    if(!all(c(a,b) %in% psi$sample)) {
        rlang::abort("Some samples in `a` and `b` are not found in `psi`.")
    }

    rlang::inform("Performing differential analysis")
    psi %>%
        dplyr::filter(sample %in% c(a,b)) %>%
        dplyr::mutate(group = ifelse(sample %in% a, "A", "B")) %>%
        dplyr::group_by(gene_id, sample, pos) %>%
        dplyr::mutate(totalcount = sum(count)) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(seqnames, start, end, strand, pos, gene_name, group) %>%
        dplyr::summarise(count = sum(count), totalcount = sum(totalcount), mean_PSI = mean(PSI)) %>%
        tidyr::unite(val, c("count", "totalcount", "mean_PSI"), sep = "_") %>%
        tidyr::spread(group, val) %>%
        dplyr::ungroup() %>%
        tidyr::separate(A, c("A.count", "A.totalcount", "A.PSI"), convert = T, sep = "_") %>%
        tidyr::separate(B, c("B.count", "B.totalcount", "B.PSI"), convert = T, "_") %>%
        dplyr::mutate(deltaPSI = B.PSI - A.PSI) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(fisher.pvalue = fisher.test(rbind(c(A.count, A.totalcount), c(B.count, B.totalcount)))$p.value) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(fisher.FDR = p.adjust(fisher.pvalue, method = "fdr")) %>%
        dplyr::select(seqnames:gene_name, A.PSI, B.PSI, deltaPSI, fisher.pvalue, fisher.FDR)
}


.bamunlist <- function (x) {
    x1 <- x[[1L]]
    if (is.factor(x1)) {
        structure(unlist(x), class = "factor", levels = levels(x1))
    } else {
            do.call(c, x)
    }
}

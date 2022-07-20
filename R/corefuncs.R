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
#' @importFrom stats end na.omit
#' @importFrom methods is
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
        dplyr::filter(afl.start==x.start | afl.end==x.end) %>%
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
#' @param colData Dataframe containing sample info
#'
#' @return
#' Data-frame containing read count and PSI for each AFL
#'
#'
.countAFL <- function(x, colData, dir, pairedEnd, n, strandSpecific) {

    count <- width <- gene_id <- type <- norm_count <- totalnormcount <- NULL
    strand <-  start <-  seqnames <-  coding <-  PSI <- NULL

    # Checks
    # catch missing args

    # retrieve list of bam files
    db <- x
    dir <- stringr::str_remove(dir, "/$")
    bams <- paste0(rownames(colData),".bam")
    if(!all(bams %in% list.files(dir))){
        rlang::abort(sprintf("BAM files not found at %s directory", dir))
    }

    dbSAF <- db %>%
        as.data.frame() %>%
        dplyr::select(GeneID=effectivecoord, Chr=seqnames, Start=start, End=end, Strand=strand)

    print(file.path(dir, bams))


    subread_out <- suppressMessages(Rsubread::featureCounts(file.path(dir, bams),
                                                            annot.ext = dbSAF,
                                                            isPairedEnd = pairedEnd,
                                                            nthreads = n,
                                                            strandSpecific = strandSpecific))
    counts <-  subread_out$counts[db$effectivecoord,]

    # normalize counts by exon length (kb)
    counts <- counts/(lengths(db)/1000)

    # get sum of all normalized counts by event
    matsum <- rowsum(counts, group = paste0(db$gene_id,"-",db$type))
    matsumfull <- matsum[paste0(db$gene_id,"-",db$type),]
    rownames(matsumfull) <- db$effectivecoord

    # rename columns
    colnames(counts) <- colnames(matsumfull) <- rownames(colData)


    out <- list(inc = counts, total = matsumfull)

    return(out)
}



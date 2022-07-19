
#' Quantify PSI metrics
#'
#' @param db rowRanges slot from Quaffle object
#' @param out counts slot from Quaffle object
#'
#' @return Matrix containing PSI values
.quantAFL <- function(db, out){

    width <- gene_id <- type <- norm_count <- totalnormcount <- NULL
    strand <- start <- coord <- PSI <- NULL

    outdf <- as.data.frame(out)
    outdf$coord <- rownames(out)

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
        dplyr::select(-norm_count) %>%
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

    genenorm <- out.comb %>%
        dplyr::select(gene_id, sample, totalnormcount) %>%
        dplyr::distinct(gene_id, sample, .keep_all = T) %>%
        tidyr::spread(sample, totalnormcount) %>%
        tibble::column_to_rownames("gene_id") %>%
        as.matrix()


    return(list(psi[db$coord,], genenorm))
}

.diff <- function(x, colData, comparison, min_samples, adjust.methods){

    id <- psilist <- B.PSI <- A.PSI <- p_val <- deltaPSI <- adj_p_val <- NULL
    grouping <- comparison[[3]]
    A <- comparison[[1]]
    B <- comparison[[2]]

    state <- x@expstate
    state <-  state %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene_id") %>%
        tidyr::gather("samples","expstate", -gene_id)
    gene.considered <- colData %>%
        dplyr::left_join(state, by = "samples") %>%
        dplyr::group_by(gene_id, !!!rlang::syms(grouping)) %>%
        #dplyr::filter(sum(expstate == "OK") >= min_samples) %>%
        dplyr::summarise(numok = sum(expstate == "OK")) %>%
        tidyr::spread(grouping, numok) %>%
        dplyr::filter((!!!rlang::syms(A)) >= min_samples & (!!!rlang::syms(B)) >= min_samples) %>%
        dplyr::pull(gene_id)

    afl.considered <- x@rowRanges[x@rowRanges$gene_id %in% gene.considered]



    # state <- state[rownames(x@c),]
    # state <- state %>%
    #     as.data.frame() %>%
    #     #dplyr::mutate(id = rownames(x)) %>%
    #     tibble::rownames_to_column("id") %>%
    #     tidyr::gather("samples","state", -id)

    psi <- .posterior(x@psi) %>%
        as.data.frame() %>%
        dplyr::mutate(id = rownames(x@psi)) %>%
        #tibble::rownames_to_column("id") %>%
        dplyr::filter(id %in% afl.considered$coord) %>%
        tidyr::gather("samples","psi", -id)





    colData %>%
        dplyr::left_join(psi, by = "samples") %>%
        # dplyr::left_join(state, by = c("samples", "id")) %>%
        dplyr::group_by(id, !!!rlang::syms(grouping)) %>%
        # dplyr::filter(sum(state == "OK") >= min_samples) %>%
        # dplyr::ungroup() %>%
        # dplyr::group_by(id, !!!rlang::syms(grouping)) %>%
        dplyr::summarise(psilist = list(psi)) %>%
        tidyr::spread(grouping, psilist) %>%
        # dplyr::filter(!S4Vectors::isEmpty(!!!rlang::syms(A)) & !S4Vectors::isEmpty(!!!rlang::syms(B))) %>%
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

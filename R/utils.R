psi <- function(x) x@assays@data$inc/x@assays@data$total

.getcoords <- function(ranges, ...){
    if(missing(...)){
        return(ranges$effectivecoord)
    } else {
        ranges %>%
            as.data.frame() %>%
            dplyr::mutate(name = gene_name, gene = gene_id, coord2 = effectivecoord) %>%
            dplyr::select(gene, name, coord2, effectivecoord) %>%
            tidyr::gather("type", "feature", -effectivecoord) %>%
            dplyr::filter(feature %in% c(...)) %>%
            dplyr::pull(effectivecoord) %>%
            na.omit() %>%
            unique()
    }
}


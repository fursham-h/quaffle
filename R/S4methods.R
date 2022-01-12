#' Quaffle Object
#'
#' @slot counts Matrix containing read counts
#' @slot rowRanges GRanges object of AFL ranges
#' @slot colData DataFrame containing sample info
#' @slot psi Matrix containing PSI values
#' @slot state Matrix containing state of counts
#'
#' @return Quaffle object
#' @export
#'
#' @examples
setClass("QuaffleObject",
         slots = c(
             counts = "matrix",
             rowRanges = "data.frame",
             colData = "data.frame",
             psi = "matrix",
             state="matrix"
         )
)



## Generics
setGeneric("head", function(x, ...) standardGeneric("head"))
setGeneric("build", function(x, ...) standardGeneric("build"))
setGeneric("count", function(x, ...) standardGeneric("count"))
setGeneric("quant", function(x, ...) standardGeneric("quant"))
setGeneric("colData", function(x) standardGeneric("colData"))
setGeneric("colData<-", function(x, value) standardGeneric("colData<-"))
setGeneric("counts", function(x) standardGeneric("counts"))
setGeneric("psi", function(x) standardGeneric("psi"))
setGeneric("diff", function(x, ...) standardGeneric("diff"))
setGeneric("ranges", function(x, ...) standardGeneric("ranges"))



setMethod("show",
          "QuaffleObject",
          function(object) {
              cat(sprintf("Quaffle object with %s samples and %s features\n", nrow(object@colData), length(object@rowRanges)))
              samples <- rownames(object@colData)
              if(length(samples) > 2){
                  samples <- c(samples[1], "...", samples[ length(samples)])
              }
              features <- object@rowRanges$coord
              if(length(features) > 2){
                  features <- c(features[1], "...", features[ length(features)])
              }
              cat(sprintf("samples (%s): %s\n", nrow(object@colData), paste0(samples, collapse = " ")))
              cat(sprintf("features (%s): %s\n", length(object@rowRanges), paste0(features, collapse = " ")))
              #cat("Preview of PSI values:\n")
              #print(object@psi)
          }
)

#setMethod("show", "QuaffleObject", function(object) object@psi)


setMethod("build", "QuaffleObject", function(x, gtf) {
    x@rowRanges <- .buildAFL(gtf)
    x
})


setMethod("count", "QuaffleObject", function(x, dir) {
    x@counts <- .countAFL(x, dir)
    x
})

setMethod("quant", "QuaffleObject", function(x, min_read) {
    x@psi <- .quantAFL(x@rowRanges, x@counts)
    x@state <- ifelse(x@counts < min_read, "LOW", "OK")
    x
})


setMethod("colData", "QuaffleObject", function(x) x@colData)
setMethod("colData<-", "QuaffleObject", function(x, value) {
    x@colData <- value
    x
})

setMethod("counts", "QuaffleObject", function(x) x@counts)
setMethod("psi", "QuaffleObject", function(x) x@psi)
setMethod("ranges", "QuaffleObject", function(x) x@rowRanges)




setMethod("diff", "QuaffleObject", function(x, contrast, min_samples = 2, methods = "fdr") {

    A.PSI <- NULL
    coldat <- x@colData
    coldat <- coldat[coldat[[contrast[3]]] %in% contrast[1:2],]
    diff <- .diff(x@psi, coldat, x@state, contrast, min_samples, methods)

    ranges <- as.data.frame(x@rowRanges)
    ranges <- ranges[,c("coord", "effectivecoord", "gene_id",
                     "gene_name", "strand", "type", "coding")]
    ranges %>% dplyr::left_join(diff, by = "coord") %>%
        dplyr::filter(!is.na(A.PSI))


})




































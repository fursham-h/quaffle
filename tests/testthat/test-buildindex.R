# test_that("Index build correctly", {
#     out <- .buildAFL(system.file("extdata", "wtap.gtf", package = "quaffle"))
#     expect_equal(BiocGenerics::start(out), c(12992309, 12992077, 12972685, 12966796))
#     expect_equal(out$type, c("AF", "AF"  , "AL", "AL" ))
#     expect_equal(out$coding, c("FALSE", "FALSE"  , "TRUE", "TRUE"))
# })

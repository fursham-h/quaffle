test_that("Index build correctly", {
    out <- buildAFLindex(system.file("extdata", "wtap.gtf", package = "AFLanalyze"))
    expect_equal(start(out), c(12994063, 12992623, 12992309, 12992077, 12985916,
                               12983457, 12981751, 12980855, 12972685, 12966796))
    expect_equal(out$pos, c("start", "end"  , "start", "start", "start",
                            "start", "end"  , "end"  , "end"  , "end"  ))
})

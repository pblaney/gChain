library(gChain)

library(testthat)



test_that('dedup() works', {

    expect_equal(dedup(c(rep(2, 10.5), rep(3, 20)))[30], "3.20")

})

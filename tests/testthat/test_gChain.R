library(gChain)

library(testthat)

library(gUtils)

















### utility functions


test_that('testing vec2ir() works', {

    expect_equal(width(vec2ir(c(1, 2, 2, 2))[1]), 2)
    expect_equal(width(vec2ir(c(1, 2, 2, 2))[2]), 1)
    expect_equal(width(vec2ir(c(1, 2, 2, 2))[3]), 1)

})



test_that('testing ir2vec() works', {

	gr = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
    expect_equal(length(ir2vec(gr)), 10)
    expect_equal(ir2vec(gr)[1], 3)
    expect_equal(ir2vec(gr)[10], 16)

})


test_that('testing squeeze() works', {

    gr = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
    expect_error(squeeze(gr))
    ir1 = IRanges(c(3,7,13), c(5,9,16))


})




test_that('seqinfo2gr() works', {

	gr = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
    expect_equal(width(seqinfo2gr(gr)), 25)
    expect_equal(length(seqinfo2gr(si)), 25)

})


## dedup 

test_that('dedup() works', {

    expect_equal(dedup(c(rep(2, 10.5), rep(3, 20)))[30], "3.20")

})

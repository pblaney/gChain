
library(gChain)
library(testthat)
library(gUtils)


Sys.setenv(DEFAULT_BSGENOME = "BSgenome.Hsapiens.UCSC.hg19::Hsapiens")


gr = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
gr2 = GRanges(1, IRanges(c(1,9), c(6,14)), strand=c('+','-'), seqinfo=Seqinfo("1", 25), field=c(1,2))
dt = data.table(seqnames=1, start=c(2,5,10), end=c(3,8,15))


test_that('testing gChain() works', {

    expect_true(is(gChain(), 'gChain'))
    expect_equal(length(gChain()), 1)
    ### dim()
    expect_equal(dim(gChain())[1], 0) 
    expect_equal(dim(gChain())[2], 0)
    ### links()
    expect_equal(width(links(gChain())$x), 0)
    expect_equal(width(links(gChain())$y), 0)
    ### values()
    expect_equal(values(gChain()), data.frame())
    ## c()
    expect_equal(values(c(gChain(gr), gChain(gr2))), data.frame())
    expect_equal(values(gChain(x = si, y = si)), data.frame())
    expect_equal(values(gChain(x = GRanges(), y = NULL)), data.frame())
    expect_equal(values(gChain(x = NULL, y = GRanges())), data.frame())
    ## Error: At least one of the objects to be concatanted is not a gChain
    expect_error(c(gChain(gr), c('foo')))



})


### gMultiply BUG
## > foo = gChain(grl.unlist(grl2))
## > foobar = foo = gChain(grl.unlist(grl1))
## > foo * foobar
## Error in cbind(1:length(gr), starts)[, 2] : subscript out of bounds
## > gMultiply(foo, foobar)
## Error in cbind(1:length(gr), starts)[, 2] : subscript out of bounds
## > 


### > expand(gChain(grl.unlist(grl2)))
### Error in gr.pad(x@.galx, cbind((1/abs.scale) * left.space, (1/abs.scale) *  : 
###   could not find function "gr.pad"



test_that('testing lift() works', {

    expect_equal(length(gChain()), 1)
    foo = gChain(grl.unlist(grl2))
    expect_equal(length(lift(foo, grl.unlist(grl2))), 0)
    expect_equal(length(lift(foo, grl.unlist(grl2), split.grl = TRUE)), 0)
    ## if (!(format %in% c('GRanges', 'df', 'df.all', 'matrix', 'GRangesList', 'trackData'))){
    expect_error(lift(foo, c(1))) ## Error in .local(.Object, x, ...) : Error: x must be Granges object
    expect_equal(length(lift(foo, grl2)), 0)
    expect_equal(length(lift(foo, IRanges(c(3,7,13), c(5,9,16)))), 0)

})


## cn

test_that('testing cn() works', {

    expect_equal(length(cn(gChain(grl2))), 0)
    test = grl2
    names(test) = NULL

})



## spChain 

test_that('testing spChain() works', {

    expect_equal(length(cn(spChain(grl2))), 0)
    test = grl2
    names(test) = NULL
    expect_equal(length(spChain(test)), 1)

})




test_that('testing spChain() works', {

    expect_equal(length(cn(gChain(grl2))), 0)

})




test_that('testing gSubset() works', {

    expect_error( gSubset(c('foo')))
    expect_equal(length(gSubset(gChain())), 1)

})




## txChain

test_that('testing duplicate() works', {

    expect_true(is(txChain(example_dnase), 'gChain'))
    
})




### duplicate

test_that('testing duplicate() works', {

    expect_true(is(duplicate(grl.unlist(grl2)), 'gChain'))
    
})





### copy

test_that('testing copy() works', {

    expect_true(is(copy(grl.unlist(grl2)), 'gChain'))
    expect_equal(length(copy(grl.unlist(grl2))), 1)
    
})




### delete

test_that('testing delete() works', {

    expect_true(is(delete(grl.unlist(grl2)), 'gChain'))
    expect_equal(length(delete(grl.unlist(grl2))), 1)
    
})



### invert

test_that('testing invert() works', {
    
    expect_true(is(invert(grl.unlist(grl2)), 'gChain'))
    expect_equal(length(invert(grl.unlist(grl2))), 1)

})



test_that('testing permute() works', {

    expect_error(permute(grl2), NA) ## check works
    expect_equal(length(permute(grl2)), 1)

})


## > rearrange(event = gr2)
## Error in `[<-`(`*tmp*`, ix.n, value = tmp[ix.n]) : 
##   subscript out of bounds

test_that('testing rearrange() works', {
    
    ## Error in gr.flipstrand(event[-1]) : Error: GRanges input only
    expect_error(rearrange(grl1))

})





### bfb()
### > bfb('chr2')
### Cycle 1 chrom widthNA MB
### 
### Error: subscript contains invalid names
### 

test_that('testing bfb() works', {

    expect_error(bfb('1'))

})





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


## > gCat(foo1, foo2)
## Error in c(x, list(...)) : 
##   Error: At least one of the objects to be concatanted is not a gChain


test_that('testing gCat() works', {

    expect_error(gCat(foo1, foo2)) ### shouldn't be an error

})



test_that('testing gUnique() works', {

    expect_true(is(gUnique(gChain()), 'gChain'))

})




test_that('testing squeeze() works', {

    gr = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
    expect_error(squeeze(gr))
    ir1 = IRanges(c(3,7,13), c(5,9,16))
    expect_equal(width(squeeze(ir1))[1], 3)
    expect_equal(width(squeeze(ir1))[3], 4)
    expect_equal(squeeze(IRanges()), IRanges())


})

### grl.split BUG; 
## ele = tryCatch(as.data.frame(grl)$element, error = function(e) e)
## NULL


## levapply()

test_that('levapply() works', {

    gr = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
    foo = levapply(width(gr), as.character(seqnames(gr)), function(x) if (length(x)>1) cumsum(c(0, x[1:(length(x)-1)])) else return(0))
    expect_equal(foo, c(0, 3, 6))

})




## seqinfo2gr()

test_that('seqinfo2gr() works', {

	gr = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
    expect_equal(width(seqinfo2gr(gr)), 25)
    expect_equal(length(seqinfo2gr(si)), 25)

})


## dedup 

test_that('dedup() works', {

    expect_equal(dedup(c(rep(2, 10.5), rep(3, 20)))[30], "3.20")

})



## gr.pad() 

test_that('gr.pad() works', {

    expect_equal(width(gr.pad(gr2, 10)[1]), 16)
    expect_equal(width(gr.pad(gr2, 10)[2]), 24)

})



## gr.refactor() 

test_that('gr.refactor() works', {

    expect_equal(as.character(seqnames(gr.refactor(gr2, sn = c('2', '2'))))[1], '2')
    expect_equal(as.character(seqnames(gr.refactor(gr2, sn = c('chr2', '2'))))[1], 'chr2')
    expect_equal(as.character(seqnames(gr.refactor(gr2, sn = c('chr2', 'foo'))))[2], 'foo')

})



test_that('gr.tostring() works', {

    expect_true(is(gr.tostring(gr2), 'character'))

})



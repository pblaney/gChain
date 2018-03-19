
library(gChain)
library(testthat)
library(gUtils)


Sys.setenv(DEFAULT_BSGENOME = "BSgenome.Hsapiens.UCSC.hg19::Hsapiens")


gr = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
gr2 = GRanges(1, IRanges(c(1,9), c(6,14)), strand=c('+','-'), seqinfo=Seqinfo("1", 25), field=c(1,2))
dt = data.table(seqnames=1, start=c(2,5,10), end=c(3,8,15))


i1 = IRanges(c(3,7,13), c(5,9,16))
i2 = IRanges(c(1,9, 32), c(6,14, 45))


test_that('testing gChain() works', {

    expect_true(is(gChain(), 'gChain'))
    expect_equal(length(gChain()), 1)
    expect_equal(width(gChain()$x), 0)
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
    expect_equal(values(c(gChain())), data.frame())
    ## Error: At least one of the objects to be concatanted is not a gChain
    expect_error(c(gChain(gr), c('foo')))
    ## check 'show()' works
    expect_error(show(gChain()), NA)
    ## if (length(x)>0 & length(y)>0 & length(x) != length(y)){
    expect_error(gChain(gr, gr2))
    ##  if (is(x, 'IRanges')){
    expect_equal(dim(gChain(x = i1))[1], 16)
    expect_equal(dim(gChain(x = i1))[2], 16)
    ##  if (is(y, 'IRanges')){
    expect_equal(dim(gChain(y = i2))[1], 45)
    expect_equal(dim(gChain(y = i2))[2], 45)
    ##  Error in validObject(.Object) : 
    ##    invalid class “gChain” object: 1: All intervals pairs must have a single scale
    ##  invalid class “gChain” object: 2: Scales not consistent with interval pair widths and .pad.left / .pad.right properties.
    expect_error(gChain(i2, i1))
    ## if (!is(x, 'GRanges') | !is(y, 'GRanges')){
    ## Error: inputs to gChain() must be of the type GRanges
    expect_error(gChain(grl2))
    expect_equal(dim(gChain(gr[1:2], gr2, val = 42))[1], 25)
    ## if (ncol(val)>0 & nrow(val)==1){
    expect_equal(dim(gChain(gr[1:2], gr2, val = data.frame(c(42))))[1], 25)
    ## if (is.vector(val)){
    expect_equal(dim(gChain(gr[1:2], gr2, val = c(42, 13)))[1], 25)
    ## x = NULL, y = NULL, pad.left = 0, pad.right = 0, scale = NULL, val = data.frame())
    ## if (is.null(scale)){
    expect_equal(dim((gChain(gr2, scale=1)))[1], 25)   
    ## what does scale < 0 mean?
    expect_equal(dim((gChain(gr2, scale=-1)))[1], 25)
    
})





test_that('testing gMultiply() works', {

    gc1 = gChain(gr, gr)
    gc2 = gChain(gr, gr)
    gc2 = gMultiply(gc1, gc2)
    ## expect_equal(gc1, gc2) 
    expect_equal(gMultiply(gc1, gc2), gMultiply(gc2, gc1)) 
    ## empty GRanges()
    foo1 = gChain(grl.unlist(grl1), grl.unlist(grl2)[1:500]) 
    foo2 = gChain(example_dnase[1:200])
    multiplied = gMultiply(foo1, foo2)
    expect_equal(multiplied, GRanges())

})





test_that('testing scale() works', {

    expect_equal(scale(gChain()), 1)

})





test_that('testing lift() works', {

    expect_equal(length(gChain()), 1)
    foo = gChain(grl.unlist(grl2))
    expect_equal(length(lift(foo, grl.unlist(grl2))), 502)
    expect_equal(length(lift(foo, grl.unlist(grl2), split.grl = TRUE)), 502)
    ## if (!(format %in% c('GRanges', 'df', 'df.all', 'matrix', 'GRangesList', 'trackData'))){
    expect_error(lift(foo, c(1))) ## Error in .local(.Object, x, ...) : Error: x must be Granges object
    expect_equal(length(lift(foo, grl2)), 251)
    expect_equal(length(lift(foo, IRanges(c(3,7,13), c(5,9,16)))), 0)
    ## data.frame
    expect_true(is(lift(foo, grl.unlist(grl2), format='df'), 'data.frame'))
    expect_equal(dim(lift(foo, grl.unlist(grl2), format='df'))[1], 502)
    expect_equal(dim(lift(foo, grl.unlist(grl2), format='df'))[2], 4)
    ## GRanges
    expect_equal(length(lift(foo, grl.unlist(grl2), format='GRanges')), 502)
    expect_true(is(lift(foo, grl.unlist(grl2), format='GRanges'), 'GRanges'))
    ## trackData
    expect_equal(dim(lift(foo, grl.unlist(grl2), format='trackData'))[1], 502)
    expect_equal(dim(lift(foo, grl.unlist(grl2), format='trackData'))[2], 4)
    ## matrix
    expect_equal(dim(lift(foo, grl.unlist(grl2), format='matrix'))[1], 502)
    expect_equal(dim(lift(foo, grl.unlist(grl2), format='matrix'))[2], 4)

})



test_that('testing scale()) works', {

    expect_equal(scale(gChain()), 1)

})





test_that('testing pads() works', {

    expect_equal(pads(gChain())$pad.left, 0)
    expect_equal(pads(gChain())$pad.right, 0)

})




test_that('testing genomes() works', {

    expect_equal(width(genomes(gChain())$x), 0)
    expect_equal(width(genomes(gChain())$y), 0)

})




test_that('testing "*" works', {

    ## identity multiplication
    ## setMethod("*", signature(e1 = "gChain", e2 = "GRanges"), function(e1, e2) 
    gc1 = gChain(gr,gr)
    gc2 = gChain(gr,gr)
    ## gc2 = gc1 * gc2
    ## expect_equal(gc1, gc2) 
    expect_equal((gc1 * gc2), (gc2 * gc1)) 
    ## with GRList
    ## setMethod("*", signature(e1 = "gChain", e2 = "GRangesList"), function(e1, e2)
    expect_equal(length(gc1 * grl2), 0)
    ## with GRanges
    expect_equal(length(gc1 * gr2), 3)

})





test_that('testing expand() works', {

    ### setMethod("expand", signature(x = "gChain"), function(x, space = NULL, shift.x = FALSE, shift.y = TRUE)
    expect_true(is(expand(gChain(grl.unlist(grl2))), 'gChain'))
    expect_equal(dim(expand(gChain(grl.unlist(grl2))))[1], 3095693983)
    expect_equal(dim(expand(gChain(grl.unlist(grl2))))[2], 3095693983)
    ##
    expect_equal(dim(expand(gChain(grl.unlist(grl2)[1:25], grl.unlist(grl1)[1:25])))[1], 3095693983)
    expect_equal(dim(expand(gChain(grl.unlist(grl2)[1:25], grl.unlist(grl1)[1:25])))[2], 3095693983)
    expect_equal(dim(expand(gChain(grl.unlist(grl2)[1:500], grl.unlist(grl1))))[1], 3095693983)
    expect_equal(dim(expand(gChain(grl.unlist(grl2)[1:500], grl.unlist(grl1))))[2], 3095693983)
    ## shift.x = FALSE, shift.y = FALSE
    expect_equal(dim(expand(gChain(grl.unlist(grl2)), shift.x=FALSE, shift.y=FALSE))[1], 3095693983)
    expect_equal(dim(expand(gChain(grl.unlist(grl2)), shift.x=FALSE, shift.y=FALSE))[2], 3095693983)
    ## shift.x = TRUE, shift.y = FALSE
    expect_equal(dim(expand(gChain(grl.unlist(grl2)), shift.x=TRUE, shift.y=FALSE))[1], 3095693983)
    expect_equal(dim(expand(gChain(grl.unlist(grl2)), shift.x=TRUE, shift.y=FALSE))[2], 3095693983)

})





test_that('testing values() works', {

    ## it looks like this is always an empty data.frame()
    expect_equal(values(gChain(grl.unlist(grl2))), data.frame())
    expect_equal(values(gChain(gr2, gr2)), data.frame())

})





test_that('testing t() works', {

    expect_true(is(t(gChain()), 'gChain'))

})





test_that('testing breaks() works', {

    ### BUG
    expect_error(breaks(gChain()))
    expect_equal(breaks(gChain(grl.unlist(grl2))), GRangesList())
    #### Error in gr.start(x@.galy, 2) - 1 : 
    ####   adjustment would result in ranges with negative widths ### happens because 'clip'
    ## en = pmin(en, end(x))
    ## st = pmax(st, start(x))
    ##
    expect_equal(breaks(gChain(grl.unlist(grl2)[1:500], grl.unlist(grl1))), GRangesList())
    ##
    expect_equal(breaks(gChain(gr, gr)), GRangesList())
    expect_equal(breaks(gChain(example_genes)), GRangesList())
    expect_error(breaks('foo'))

})




## cn
test_that('testing cn() works', {

    expect_equal(length(cn(gChain(grl.unlist(grl2)))), 1029)
    test = grl2
    names(test) = NULL
    expect_equal(length(cn(gChain(grl.unlist(grl2)))), 1029) 
    expect_equal(length(cn(gChain(grl.unlist(grl2)), rev=TRUE)), 1029)

})


## setMethod('[', 'gChain', function(x, i){
test_that('testing [ works', {

    expect_true(is(gChain()[1], 'gChain'))
    ## Error: subscript contains out-of-bounds indices
    expect_error(gChain()[2])

})


## spChain 
test_that('testing spChain() works', {

    expect_true(is(spChain(grl2[1:3]), 'gChain'))
    expect_equal(width(si2gr(seqinfo(spChain(grl2[1:3]))$y)[1]), 2)
    ##  if (is.null(names(grl))){
    foo = grl2
    names(foo) = NULL
    expect_true(is(spChain(foo[1:3]), 'gChain'))
    expect_equal(width(si2gr(seqinfo(spChain(foo[1:3]))$y)[1]), 2)   

})





test_that('testing gSubset() works', {

    expect_error( gSubset(c('foo')))
    expect_equal(length(gSubset(gChain())), 1)

})



## paChain()

test_that('testing paChain() works', {

    s1 = DNAString("ACTTCACCAGCTCCCTGGCGGTAAGTTGATCAAAGGAAACGCAAAGTTTTCAAG")
    s2 = DNAString("GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC")
    expect_equal(dim(paChain('ACTTCACCAGCTCCCTGGCGGTAAGTTGATCAAAGGAAACGCAAAGTTTTCAAG', 'GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC'))[1], 54)
    expect_equal(dim(paChain('ACTTCACCAGCTCCCTGGCGGTAAGTTGATCAAAGGAAACGCAAAGTTTTCAAG', 'GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC'))[2], 60)   
    expect_equal(dim(suppressWarnings(paChain(s1, s2)))[1], 54)
    expect_equal(dim(suppressWarnings(paChain(s1, s2)))[2], 60)
    ## if (inherits(seq1, 'factor')){
    factor1 = factor("ACTTCACCAGCTCCCTGGCGGTAAGTTGATCAAAGGAAACGCAAAGTTTTCAAG")
    factor2 = factor("GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC")
    expect_equal(dim(suppressWarnings(paChain(factor1, factor2)))[1], 54)
    expect_equal(dim(suppressWarnings(paChain(factor1, factor2)))[2], 60)
    ## score.thresh
    expect_equal(dim(suppressWarnings(paChain(factor1, factor2, score.thresh = 5000)))[1], 54)
    expect_equal(dim(suppressWarnings(paChain(factor1, factor2, score.thresh = 5000)))[2], 60)   
    # reference sequence
    # standardseq = AAString("MARKSLEMSIR")
    # query sequences
    # seq1 = AAString("MARKSLEMSER")
    # seq2 = AAString("MDRKSLEMSER")
    # seq3 = AAString("MDRKSAEMSER")
    # seq4 = AAString("MARKSLEMSIR")
    proteinseq1 = AAString("MDRKSAEMSER")
    proteinseq2 = AAString("MDRKSLEMSER")
    expect_equal(dim(suppressWarnings(paChain(proteinseq1, proteinseq2)))[1], 11)
    expect_equal(dim(suppressWarnings(paChain(proteinseq1, proteinseq2)))[2], 11)
    ## assume.dna=FALSE
    expect_equal(dim(paChain('MDRKSAEMSER', 'MDRKSLEMSER', assume.dna=FALSE))[1], 11)
    expect_equal(dim(paChain('MDRKSAEMSER', 'MDRKSLEMSER', assume.dna=FALSE))[2], 11)
    ## Error in .Call2("new_XString_from_CHARACTER", classname, x, start(solved_SEW),  : 
    ##   key 69 (char 'E') not in lookup table
    expect_error(dim(paChain('MDRKSAEMSER', 'MDRKSLEMSER')))    
    expect_error(paChain(AAString("MDRKSAEMSER"), 'ACCCT'))
    ## 
    fact1 = "ACTTCACCAGCTCCCTGGCGGT"
    fact2 = "GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC"
    expect_equal(dim(suppressWarnings(paChain(fact1, fact2)))[1], 22)
    expect_equal(dim(suppressWarnings(paChain(fact1, fact2)))[2], 60)
    fa1 = "A"
    fa2 = "GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC"
    expect_equal(dim(suppressWarnings(paChain(fa1, fa2)))[1], 1)
    expect_equal(dim(suppressWarnings(paChain(fa1, fa2)))[2], 60)
    f1 = "ACTTCACCAGCTCCCTGGCGGT"
    f2 = "C"
    expect_equal(dim(suppressWarnings(paChain(f1, f2)))[1], 22)
    expect_equal(dim(suppressWarnings(paChain(f1, f2)))[2], 1)
    ##  *** caught segfault ***
    ## address 0x7f7fc7800000, cause 'memory not mapped'
    ## pairwiseAlignment(seq1, seq2)
    ff1 = "ACTTCACCAGCTCCCTGGCGGT"
    ff2 = ""
    ##     Error in paChain(ff1, ff2) : 
    ##   Either seq1 and/or seq2 is of width 0. This will result in an error with pairwiseAlignment(seq1, seq2).
    expect_error(paChain(ff1, ff2))
    ##
    fact1 = "ACTTCACCAGCTCCCTGGCGGT"
    fact2 = "GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC"
    expect_equal(dim(suppressWarnings(paChain(fact1, fact2, sn1 = c('one', 'two', 'three', 'four'), sn2 = c('one', 'one', 'one'))))[1], 60)
    expect_equal(dim(suppressWarnings(paChain(fact1, fact2, sn1 = c('one', 'two', 'three', 'four'), sn2 = c('one', 'one', 'one'))))[2], 60)

})







## cgChain()

test_that('testing duplicate() works', {

    ## default 
    gChain_duplicate = duplicate(gr2)
    expect_equal(dim(gChain_duplicate)[1], 25)
    expect_equal(dim(gChain_duplicate)[2], 37)
    foo = GRanges(1, IRanges(c(1, 8), c(10, 14)), strand=c('+','+'))
    ## Error in duplicate(foo) : Error: Input ranges must be nonoverlapping
    expect_error(duplicate(foo))
    expect_error(duplicate(example_genes))    
    
})



## maChain()


## txChain

test_that('testing duplicate() works', {
    
    ## GRangesList input
    expect_equal(dim(txChain(grl2))[1], 3095693983)
    expect_equal(dim(txChain(grl2))[2], 502)
    ## GRanges input
    expect_true(is(txChain(example_dnase), 'gChain')) 
    ## val
    labels = c(rep('hey', 10000))
    expect_equal(dim(txChain(example_dnase, val=labels))[1], 3095693983)
    expect_equal(dim(txChain(example_dnase, val=labels))[2], 6070884)
    ##  Error: val data.frame must correspond to grl, i.e. have nrows = length(grl)
    expect_error(txChain(example_dnase, val = c("foo", "bar")))
    ## txchain
    txchain_labels = c(rep('foo', 251))
    expect_equal(dim(txChain(example_dnase, txchain_labels))[1], 3095693983)
    expect_equal(dim(txChain(example_dnase, txchain_labels))[2], 6070884)
    ## exonFrame.field 
    expect_true(is(txChain(grl2, exonFrame.field = 'bin'), 'gChain')) 
    expect_equal(dim(txChain(grl2, exonFrame.field = 'bin'))[1], 3095693983)
    expect_equal(dim(txChain(grl2, exonFrame.field = 'bin'))[2], 761)  
    ## translate
    expect_equal(dim(txChain(example_dnase, translate=TRUE))[1], 3095693983)
    expect_equal(dim(txChain(example_dnase, translate=TRUE))[2], 2023628)

})




### duplicate

test_that('testing duplicate() works', {

    expect_true(is(duplicate(grl.unlist(grl2)), 'gChain'))
    
})





### copy

test_that('testing copy() works', {

    expect_true(is(copy(grl.unlist(grl2)), 'gChain'))
    expect_equal(length(copy(grl.unlist(grl2))), 1)
    ## if (any(strand(from) == '*')){
    foo = example_dnase
    strand(foo) = '*'
    expect_true(is(copy(foo), 'gChain'))

    
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
    ##  if (any(strand(c.gr)=='*'))
    foo = grl2
    strand(foo) = '*'


})


## > rearrange(event = gr2)
## Error in `[<-`(`*tmp*`, ix.n, value = tmp[ix.n]) : 
##   subscript out of bounds

test_that('testing rearrange() works', {
    
    ## Error in gr.flipstrand(event[-1]) : Error: GRanges input only
    expect_error(rearrange(grl1))

})





test_that('testing bfb() works', {

    ## stop('w and wd must be positive')
    expect_error(bfb(w=-.05))
    expect_error(bfb(3), NA) ## check runs
    ### even setting set.seed(), the dim() varies
    expect_true(is(bfb(3), 'gChain'))

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
    expect_equal(length(ir2vec(i2, each=c(2, 4, 6))), 120)

})


## > gCat(foo1, foo2)
## Error in c(x, list(...)) : 
##   Error: At least one of the objects to be concatanted is not a gChain

test_that('testing gCat() works', {

    expect_true(is(gCat(), 'gChain'))
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
    irlength1 = IRanges(c(1300), c(1600))
    expect_equal(as.data.frame(squeeze(irlength1))$width, 301)



})

### grl.split BUG; 
## ele = tryCatch(as.data.frame(grl)$element, error = function(e) e)
## NULL
test_that('grl.split() works', {

    expect_error(grl.split(grl2))

})



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
    expect_equal(length(seqinfo2gr(si, strip.empty=TRUE)), 25)

})


## dedup()
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
    ## if (is.factor(sn)){
    expect_equal(as.character(seqnames(gr.refactor(gr2, sn = factor(c('2', '2')))))[1], '2')        
    expect_equal(as.character(seqnames(gr.refactor(gr2, sn = factor(c('chr2', '2')))))[1], 'chr2')
    expect_equal(as.character(seqnames(gr.refactor(gr2, sn = factor(c('chr2', 'foo')))))[2], 'foo')

})


## gr.tostring()
test_that('gr.tostring() works', {

    expect_true(is(gr.tostring(gr2), 'character'))

})




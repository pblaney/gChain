###############################################################################
## Marcin Imielinski
## Weill Cornell Medicine mai9037@med.cornell.edu
## NYGC mimielinski@nygenome.org

## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.

## You should have received a copy of the GNU Lesser General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################


########################
# class::gChain
#
# Class representing genomic coordinate transformations resulting from alignments, gene-protein mappings, and other homologous
# or syntenic relationships between sequences.
#
#
## Marcin Imielinski
## The Broad Institute of MIT and Harvard / Cancer program.
## marcin@broadinstitute.org
########################



#################################################################
#' @name gChain-class
#' @title gChain-class
#' @description
#'
#' specified by pair of GRanges objects .galx and .galy, each defined on a genome (i.e. a Seqinfo object
#' specifying seqlevels and seqlengths, contained inside a GRanges object).  These objects are the
#' same length and each pair of ranges .galx[i], .galy[i] specifies a region of synteny on each "genome".
#'
#' The widths of mapped interval pairs must be consistent with the (1) scale s and (2) .frac.left / .frac.right parameters
#' For non fractional mappings, the scale is the ratio of widths of the y and x interval pairs, which the same number
#' for every pair.   The scale must be an integer or the reciprocal of an integer.
#'
#' For scales s != 1 (eg gene to protein mappings) fractional mappings (e.g. exon to protein) may result in which
#' the first position on the y interval is mapped to by only k* < k = s (or 1/s) coordinates on the y interval
#' (for a given pair).   This occurs with exon intervals whose boundaries lie inside a codon.  
#'
#' For this we allow "padded" mappings where each interval pair is associated with a pad.left and pad.right,
#' which are offsets between 0 and k-1 (where k is the scale or 1/scale, whichever is the integer).
#' i.e. if scale > 1,
#' then k-pad.left specifies the number of coordinates on the left side of the x interval that will map to the
#' to the leftmost coordinate on the y interval.
#' then k-pad.right specifies the number of coordinates on the right side of the x interval that will map to the
#' to the rightmost coordinate on the y interval.
#' (for 0 < scale < 1, just switch x and y)
#'
#' i.e. if scale < -1,
#' then k-pad.left specifies the number of coordinates on the LEFT side of the x interval that will map to the
#' to the RIGHTmost coordinate on the y interval.
#' then k-pad.right specifies the number of coordinates on the RIGHT side of the x interval that will map to the
#' to the LEFTmost coordinate on the y interval.
#' (for -1 < scale < 0, just switch x and y)
#'
#' ie padding always refers to genome with ranges that are being contracted in the chain
#' 
#' These interval pairs won't have widths that are integer multiples, however (for scales>1) the  (width(y)+pad.left+pad.right)/width(x) 
#' will be an integer, or (for scales<1)  (width(x)+pad.left+pad.right)/width(y) will be an integer.
#' 
#' @import gUtils GenomicRanges Matrix
#################################################################
setClass('gChain', representation(.galx = 'GRanges', .galy = 'GRanges', .scale = 'numeric', .pad.left = 'integer', .pad.right = 'integer', values = 'data.frame', .n = 'numeric', .m = 'numeric'))

suppressWarnings(removeMethod('show', 'gChain')) 

setMethod('show', 'gChain', function(object){ 
    message(sprintf('gChain object with scale(s) %s mapping %s GRanges on sequence of length %s to %s GRanges on sequence of length %s.\n',
        paste(unique(object@.scale), collapse = ", "), length(object@.galx), object@.n, length(object@.galy), object@.m))
})




#' @name initialize
#' @title initialize
#' @description
#'
#' gchain::initialize
#' gchain initialize
#'
#' if only x (or y) GRanges is given will make "identity chain"
#' if x (or y) is given as a seqinfo, then will make "empty chain"
#'
#' @export
setMethod('initialize', 'gChain', function(.Object, x = NULL, y = NULL, pad.left = 0, pad.right = 0, scale = NULL, val = data.frame())
{
    
    if (length(x)>0 & length(y)>0 & length(x) != length(y)){
        stop('x and y GRanges must be of the same length')
    }
  
    if (is(x, 'Seqinfo')){
        x = GRanges(seqlengths = seqlengths(x))
    }

    if (is(y, 'Seqinfo')){
        y = GRanges(seqlengths = seqlengths(y))
    }
            
    if (is.null(x) & is.null(y)){
        x = GRanges('NA', IRanges(1, 0))
        y = GRanges('NA', IRanges(1, 0))
    } else{
        ## if x or y undefined then we just make this the "identity chain", ie filling in the other
        if (is.null(x)){
            x = y
        } else if (is.null(y)){
            y = x
        }
    }
            
    if (is(x, 'IRanges')){
        x = GRanges('NA', x)
    }

    if (is(y, 'IRanges')){
        y = GRanges('NA', y)
    }
                          
    if (any(is.na(seqlengths(x)))){
        x = gr.fix(x, drop = F)
    }
    if (any(is.na(seqlengths(y)))){
        y = gr.fix(y, drop = F)
    }

    if (!is(x, 'GRanges') | !is(y, 'GRanges')){
        stop('Error: inputs to gChain() must be of the type GRanges')   
    }

    keep = which(width(x)!=0 & width(y)!=0)
    if (any(!keep)) {
        x = x[keep]
        y = y[keep]
    }
            
    strand(x)[which(as.logical(strand(x)=='*'))] = '+'                                                                                                                  
    strand(y)[which(as.logical(strand(y)=='*'))] = '+'   
            

    if (ncol(mcols(x))){
        .Object@.galx = x[, c()]
    } else{
        .Object@.galx <- x
    }

    if (ncol(mcols(y))){
        .Object@.galy = y[, c()]
    } else{
        .Object@.galy <- y
    }
  
    .Object@.n = sum(as.numeric(seqlengths(x)), na.rm = T);
    .Object@.m = sum(as.numeric(seqlengths(y)), na.rm = T);
    .Object@.pad.left = as.integer(cbind(1:length(x), pad.left)[,2]) 
    .Object@.pad.right = as.integer(cbind(1:length(x), pad.right)[,2])

    strand.match = (strand(x)==strand(y))
    if (is.null(scale)){
        if (all(width(y)==0)){
            .Object@.scale = 1
        } else if (all((width(x)/width(y))==1)){
            .Object@.scale = c(-1, 1)[1+as.numeric(strand.match)]                  
        } else if (all((width(y)/width(x))>=1)){
            .Object@.scale = (width(y) + pad.left + pad.right)/width(x)*c(-1, 1)[1+as.numeric(strand.match)]
        } else if (all((width(x)/width(y))>=1)){
            .Object@.scale = width(y)/(width(x) + pad.left + pad.right)*c(-1, 1)[1+as.numeric(strand.match)]
        } else{
            stop('Error: Ambiguous scale specified by widths.  Please provide directly')
        }
    } else{
        .Object@.scale = scale[1]*c(-1, 1)[1+as.numeric(strand(x)==strand(y))]
    }

    if (is.null(val)){
        val = data.frame()
    }

    if (!is.data.frame(val)){
        val = as.data.frame(val)
    }

    if (is.vector(val)){
        val = data.frame(val = val)
    }
            
    if (ncol(val)>0 & nrow(val)==1){
        .Object@values = do.call('rbind', lapply(1:length(.Object@.galx), function(y) val))
    } else{
        .Object@values = val;
    }

    validObject(.Object)            
    return(.Object)          
})


setValidity('gChain', function(object){

    problems = c();

    if (length(object@.galx) == length(object@.galy)){
        if (any(width(object@.galx)>0)){
            abs.scale = unique(abs(object@.scale))

            if (length(abs.scale)>1){
                problems = c(problems, 'All intervals pairs must have a single scale');                
            }
                    
            if (!((abs.scale %% 1)==0 | ((1/abs.scale) %% 1)==0)){
                problems = c(problems, 'Interval pairs must have a scale that is either an integer or the reciprocal of an integer.');
            }

            if (abs.scale > 1){                        
                if (!all((width(object@.galy) + object@.pad.left + object@.pad.right) == width(object@.galx)*abs.scale)){
                    problems = c(problems, 'Scales not consistent with interval pair widths and .pad.left / .pad.right properties.');
                }

                if (any(object@.pad.left>=abs.scale) | any(object@.pad.right>=abs.scale)){
                    problems = c(problems, '.pad.left and .pad.right must be strictly smaller than scale for scales greater than 1')
                }
            } else{
                if (!all(((width(object@.galx) + object@.pad.left + object@.pad.right)*abs.scale) == width(object@.galy))){
                    problems = c(problems, 'Scales not consistent with interval pair widths and .pad.left / .pad.right properties.');
                }

                if (any(object@.pad.left>=(1/abs.scale)) | any(object@.pad.right>=(1/abs.scale))){
                    problems = c(problems, '.pad.left and .pad.right must be strictly smaller than  1/scale for scales less than 1');      
                }
            }
        }           
    } else{
        problems = c(problems, 'Length of interval pairs do not match')
    }

    if (nrow(object@values)>0){
        if (nrow(object@values) != length(object@.galx)){
            problems = c(problems, 'Length of values incompatible with gChain intervals')
        }
    }
                
    if (length(problems)==0){
        TRUE
    } else{
        problems
    }
})







#' @name gMultiply
#' @name gMultiply
#' @export
setGeneric('gMultiply', function(e1, e2, pintersect=NA) standardGeneric('gMultiply'))
setMethod("gMultiply", signature(e1 = "gChain", e2 = "gChain"), function(e1, e2, pintersect=NA){

    ##if (!.identical.seqinfo(seqinfo(e2)$y, seqinfo(e1)$x))
    ##  warning('Genomes of range of e2 of domain of e1 are not identical');

    # image of e2 y intervals under e1, (e2y.image is in genome C)
    e2y.image <- lift(e1, e2@.galy, pintersect=pintersect) #, mc.cores=mc.cores) #max.chunk=1e7, mc.cores=mc.cores))
    if (length(e2y.image)==0){
        return(GRanges())
    } else{
        # e2 y intervals trimmed to align with e2y.image, (e2y.preimage is in genome B)
        e2y.preimage = gr.trim(e2@.galy[values(e2y.image)$query.id], values(e2y.image)$query.start, values(e2y.image)$query.end);

        # keep track of e2 link.id (ie the query id from the first lift)
        values(e2y.preimage)$e2.link.id = values(e2y.image)$query.id

        # trimmed e2 y intervals lifted backward through e2 (e2x.preimage is in genome A)
        e2x.preimage <- lift(t(e2), e2y.preimage, pintersect=pintersect) #, mc.cores=mc.cores)   
 
        # only keep mappings that have e2.link.id = link.id
        # (avoid redundant mappings arising from back and forth lift)
        e2x.preimage = e2x.preimage[values(e2x.preimage)$e2.link.id == values(e2x.preimage)$link.id]
    
        # replicate e2y image and preimage along backlifted hits
        e2y.preimage = e2y.preimage[values(e2x.preimage)$query.id]
        e2y.image = e2y.image[values(e2x.preimage)$query.id]
  
        new.scale = e1@.scale[e2y.image$link.id]*e2@.scale[e2x.preimage$link.id]

        ## save link ids for downstream processing
        e1.link.ids = values(e2y.image)$link.id
        e2.link.ids = values(e2x.preimage)$link.id
  
        pad.left = 0;
        pad.right = 0;
  
        if (length(new.scale)>0){    
            if (all(abs(new.scale)>1)){
                ## to determine padding on right and left 
                ## need to traverse both chains forward, lifting just the starts and end points
                ## of each e1x.preimage interval and noting how many e2y.image coordinates it maps to
                ## if this is equal to new.scale on each side then no pad necessary,
                ## otherwise we need padding

                # keeping track of link.ids is crucial to make sure we don't do redundant traversals
                # if chains are very tangled / promiscuous
                values(e2x.preimage)$e1.link.id = e1.link.ids;
                values(e2x.preimage)$e2.link.id = e2.link.ids;
        
                e2y.starts = lift(e2, gr.trim(e2x.preimage, 1))
                e2y.starts = e2y.starts[values(e2y.starts)$link.id == values(e2y.starts)$e2.link.id]
                e2y.starts = e2y.starts[order(values(e2y.starts)$query.id)]
            
                e2y.ends = lift(e2,  gr.trim(e2x.preimage, start = width(e2x.preimage)))
                e2y.ends = e2y.ends[values(e2y.ends)$link.id == values(e2y.ends)$e2.link.id]
                e2y.ends = e2y.ends[order(values(e2y.ends)$query.id)]

                ## then e1
                e2y.starts = lift(e1, e2y.starts)
                e2y.starts = e2y.starts[values(e2y.starts)$link.id == values(e2y.starts)$e1.link.id]
                e2y.starts = e2y.starts[order(values(e2y.starts)$query.id)]
        
                e2y.ends = lift(e1,  e2y.ends)
                e2y.ends = e2y.ends[values(e2y.ends)$link.id == values(e2y.ends)$e1.link.id]
                e2y.ends = e2y.ends[order(values(e2y.ends)$query.id)]
        
                pad.right = abs(new.scale)[1] - width(pintersect(e2y.ends, e2y.image))
                pad.left = abs(new.scale)[1] - width(pintersect(e2y.starts, e2y.image))

                ## flip pad right and pad left for neg sccale
                if (any(neg.map <- new.scale<0)){
                    tmp = pad.right[neg.map]
                    pad.right[neg.map] = pad.left[neg.map]
                    pad.left[neg.map] = tmp[neg.map]
                }
            
                # one codon edge case - padding is meaningless, but to preserve "scale"
                # we want to make padding + width consistent while making pad.left 0 for new.scale<0 links
                # and pad.right 0 for new.scale > 0 one.codon links
                one.codon = width(e2x.preimage)==1
                pad.left[one.codon & new.scale>0] = 0
                pad.right[one.codon & new.scale<0] = 0
            } else if (all(abs(new.scale)<1)) {
                ## now do the same in reverse if the scale is changing in the other direction
                values(e2y.image)$e1.link.id = e1.link.ids
                values(e2y.image)$e2.link.id = e2.link.ids
        
                e1x.starts = lift(t(e1), gr.trim(e2y.image, 1))
                e1x.starts = e1x.starts[values(e1x.starts)$link.id == values(e1x.starts)$e1.link.id]
                e1x.starts = e1x.starts[order(values(e1x.starts)$query.id)]
    
                e1x.starts = lift(t(e2), e1x.starts)
                e1x.starts = e1x.starts[values(e1x.starts)$link.id == values(e1x.starts)$e2.link.id]
                e1x.starts = e1x.starts[order(values(e1x.starts)$query.id)]
        
                e1x.ends = lift(t(e1), gr.trim(e2y.image, start = width(e2y.image)))
                e1x.ends = e1x.ends[values(e1x.ends)$link.id == values(e1x.ends)$e1.link.id]
                e1x.ends = e1x.ends[order(values(e1x.ends)$query.id)]
        
                e1x.ends = lift(t(e2), e1x.ends)
                e1x.ends = e1x.ends[values(e1x.ends)$link.id == values(e1x.ends)$e2.link.id]
                e1x.ends = e1x.ends[order(values(e1x.ends)$query.id)]
            
                pad.left = (1/abs(new.scale)[1])-width(pintersect(e1x.starts, e2x.preimage))
                pad.right = (1/abs(new.scale)[1])-width(pintersect(e1x.ends, e2x.preimage))

                ## flip pad right and pad left for neg sccale
                if (any(neg.map <- new.scale<0)){    
                    tmp = pad.right[neg.map]
                    pad.right[neg.map] = pad.left[neg.map]
                    pad.left[neg.map] = tmp[neg.map]
                }
        
                ## one codon edge case - padding can be on either side, we choose right
                one.codon = width(e2y.image)==1
                pad.left[one.codon & new.scale>0] = 0
                pad.right[one.codon & new.scale<0] = 0
            }
        }
  
        val.e1 = values(e1)[e2y.image$link.id, ,drop = FALSE]
        val.e2 = values(e2)[e2x.preimage$link.id, ,drop = FALSE]

        if (ncol(val.e1)>0 & ncol(val.e2)>0 && FALSE) {
            val = cbind(val.e1[, setdiff(names(val.e1), names(val.e2)), drop = FALSE], val.e2[, setdiff(names(val.e2), names(val.e1)), drop = FALSE])
            shared.names = intersect(names(val.e2), names(val.e1))
            if (length(shared.names)>0){
                val[, shared.names] = val.e2[,shared.names] | val.e1[,shared.names]
            }
        } else if (ncol(val.e1)>0){
            val = val.e1
        } else{
            val = val.e2
        }

        return(gChain(e2x.preimage, e2y.image, pad.left = pad.left, pad.right = pad.right, val= val))
    
    }

})




###
#' @name gMultiply
#' @title gMultiply
#' @description
#'
#' gChain multiply
#' gchain::multiply
#'
#' if e2 lifts from genome A to genome B
#' and e1 lifts from genome B to genome C
#' then e1 * e2 = lifts from genome A to genome C
#  @export
###
setMethod("*", signature(e1 = "gChain", e2 = "gChain"), function(e1, e2) {
    return(gMultiply(e1, e2, pintersect = NA))
})


setMethod("*", signature(e1 = "gChain", e2 = "GRanges"), function(e1, e2) {
    return(lift(e1, e2))
})

setMethod("*", signature(e1 = "gChain", e2 = "GRangesList"), function(e1, e2) {
    return(lift(e1, e2))
})


setGeneric('lift', function(.Object, x, ...) standardGeneric('lift')) 




setGeneric('lift', function(.Object, x, ...) standardGeneric('lift')) 

#' @name lift
#' @title lift
#' @description
#'
#' Takes GRanges input and returns GRanges output
#' Output GRanges values() fields are tagged with the following values to allow mapping to input
#' (1) query.id - index of query matching this GRanges
#' (2) query.start - start of this GRanges in the query.id query
#' (3) query.end - end of this GRanges in the query.id query
#' (4) link.id - index of link in gChain used to make this mapping
#'
#' Output GRanges will retain all other values features of the query gRanges that yielded it.
#'
#' GRangesList can be also provided as input x, in which case the output will be a GRL or list of data frames
#' one for each input
#'
#' If x is a gTrack object, then output will be a gTrack object.  gTrack inputs should be of length 1.
#'
#' if split.grl = TRUE, grl outputs are split via grl.split according to (mapped) seqname and strand
#'
#' Remaining args passed on to gr.findoverlaps
#'
#' @param x GRanges to lift through this chain
#' @param split.grl flag whether or not to split output into GRangesList for GRangesList input
#' @param format INFO INFO
#' @author Marcin Imielinski
#' @export
setMethod('lift', signature('gChain'), function(.Object, x, format = 'GRanges', split.grl = FALSE, pintersect = NA, by = NULL, verbose=TRUE, ...){

            if (!(format %in% c('GRanges', 'df', 'df.all', 'matrix', 'GRangesList', 'gTrack', 'data.frame'))){
                stop('Output format can only be "GRanges", "GRangesList", "gTrack", "data.frame", df", "df.all",  or "matrix"')
            }

            if (is(x, 'gTrack')){
                if (length(x)>1){
                  return(do.call('c', lapply(1:length(x), function(i) lift(.Object, x[i]))))
                }
                
                x.track = x[1];
                x = x.track@data[[1]]
                format = class(x);
            } else{
                x.track = NULL;
            }
            
            if (is(x, 'GRangesList')){
                input.grl = TRUE
                grl.names = names(x);
                if (is.null(grl.names)){
                    grl.names = as.character(1:length(x))
                }
                
                grl.val = values(x);
                rownames(grl.val) = grl.names
                

                tmp.df = tryCatch(as.data.frame(x), error = function(e) e)

                if (!inherits(tmp.df, 'error')){
                    #list.id = tmp.df$element
                    list.id = tmp.df$group
                    gr.name = rownames(tmp.df);
                } else {
                    ## gr names are screwy so do some gymnastics
                    if (!is.null(names(x))){
                        x.name = names(x)
                    } else{
                        x.name = 1:length(x)                    
                    }
                    list.id = as.character(Rle(x.name, sapply(x, length)))
                    tmp.x = x;
                    names(tmp.x) = NULL;
                    tmp.x = unlist(tmp.x)
                    gr.name = names(tmp.x);
                }
                
                x = unlist(x)
                values(x)$list.id = list.id;
                names(x) = gr.name;
                
                format = 'GRanges';
            } else{
                input.grl = FALSE
            }
            
            if (is(x, 'IRanges')){
                x = GRanges('NA', x)
            }

            if (!inherits(x, 'GRanges')){
                stop('Error: x must be Granges object')              
            }

            if (length(.Object@.galx)>0){

                pval <- length(seqlevels(x)) > 50 && length(seqlevels(.Object@.galx)) > 50
                if (verbose){
                    print(paste('psmart is', pval))
                }
                ## if (!is.na(pval) && is.na(pintersect))
                ##   hits <- gr.findoverlaps(x, .Object@.galx, verbose = verbose, ...) # pairs of matches
                ## else

                if (is.null(by)){
                    system.time(hits <- gr.findoverlaps(x, .Object@.galx, verbose = verbose,  ...))
                } else{
                    tmpx = .Object@.galx
                    values(tmpx) = values(.Object)
                    system.time(hits <- gr.findoverlaps(x, tmpx, verbose = verbose, by = by,   ...))
                }

                s.overlap = .Object@.scale[values(hits)$subject.id]
                s.abs = abs(s.overlap)
                neg.map = s.overlap<0
                qd.overlap = ranges(hits)
                qr.hits = .Object@.galy[values(hits)$subject.id];
                
                link.starts = start(.Object@.galx)[values(hits)$subject.id]
                link.ends = end(.Object@.galx)[values(hits)$subject.id]
                starts = ends = rep(NA, length(qr.hits))
                
                ## lift onto range
                if (unique(abs(.Object@.scale))<1){
                    pad.left.x = .Object@.pad.left[values(hits)$subject.id]
                    pad.right.x = .Object@.pad.right[values(hits)$subject.id]                    
                    
                    starts[neg.map] = start(qr.hits)[neg.map] + ceiling((link.ends[neg.map] + pad.right.x[neg.map] + 1 - end(qd.overlap)[neg.map])*s.abs[neg.map]) - 1
                    ends[neg.map] = start(qr.hits)[neg.map] + ceiling((link.ends[neg.map] + pad.right.x[neg.map]  + 1 - start(qd.overlap)[neg.map])*s.abs[neg.map]) - 1

                    starts[!neg.map] = start(qr.hits)[!neg.map] + ceiling((start(qd.overlap)[!neg.map] - link.starts[!neg.map] + 1 + pad.left.x[!neg.map])*s.abs[!neg.map]) - 1 
                    ends[!neg.map] = start(qr.hits)[!neg.map] + ceiling((end(qd.overlap)[!neg.map] - link.starts[!neg.map] + 1 + pad.left.x[!neg.map])*s.abs[!neg.map]) - 1 
                } else {
                    pad.left.y = .Object@.pad.left[values(hits)$subject.id]
                    pad.right.y = .Object@.pad.right[values(hits)$subject.id]
                    
                    shift1 = (start(qd.overlap) - link.starts)*abs(s.overlap)
                    shift2 = ((end(qd.overlap) - link.starts + 1)*(abs(s.overlap)))-1

                    starts[neg.map] = end(qr.hits)[neg.map] - (shift2[neg.map] - pad.right.y[neg.map])
                    ends[neg.map] = end(qr.hits)[neg.map] - (shift1[neg.map] - pad.right.y[neg.map])
                    starts[!neg.map] = start(qr.hits)[!neg.map] + (shift1[!neg.map] - pad.left.y[!neg.map])
                    ends[!neg.map] = start(qr.hits)[!neg.map] + (shift2[!neg.map] - pad.left.y[!neg.map])
                }
                out = GRanges(seqnames(qr.hits), IRanges(starts, ends), strand = strand(qr.hits), seqlengths = seqlengths(qr.hits));
                
                ## propagate strand flips depending on sign of s.overlap for link pair            
                has.strand.x = as.logical(as.character(strand(x)) %in% c('+', '-'))[values(hits)$query.id]
                flip=  has.strand.x  & s.overlap<0
                if (any(flip)){
                    strand(out)[flip] = c('-', '+')[1 + as.logical(strand(x)=='-')][values(hits)$query.id][flip]
                }
                if (any(!flip)){
                    strand(out)[!flip] = strand(x)[values(hits)$query.id][!flip]
                }

                ## if y links have unspecified strand then do not propagate strand information 
                has.strand.y = as.logical(as.character(strand(.Object@.galy)) %in% c('+', '-'))[values(hits)$subject.id]
                if (any(!has.strand.y)){
                    strand(out)[!has.strand.y] = '*';
                }
                
                ## propagate query GRanges values
                values(out) = values(x)[values(hits)$query.id, , drop = FALSE]
                
                ## save query indices and coordinates
                values(out)$query.id = values(hits)$query.id
                values(out)$query.start = start(qd.overlap)-start(x)[values(hits)$query.id] + 1
                values(out)$query.end = end(qd.overlap) - start(x)[values(hits)$query.id] + 1            
                values(out)$link.id = values(hits)$subject.id

                if (length(out)>0){
                    out = out[order(values(out)$query.id, values(out)$query.start)]
                }

            } else{
                out = GRanges(seqlengths = seqlengths(.Object@.galy))
            }
            ## we will expand the output 
            if (format != 'GRanges'){
                query.widths = width(x[values(out)$query.id])
                query.length = sum(query.widths)
                query.offsets = c(0, cumsum(query.widths[1:(length(query.widths)-1)]))[1:length(query.widths)]
                link.scales = .Object@.scale[values(out)$link.id]
                x.expand = pmax(abs(link.scales), 1)
                
                ## skeleton df to catch all coordinates of all hits, expanding for links with abs(scales) > 1
                df.out = data.frame(query.coord = 1:query.length, query.id = as.vector(Rle(values(out)$query.id, query.widths)),
                  chr = NA, pos = NA, stringsAsFactors = F)
                
                ## identify (flattened) x indices that map to at least one y index
                mapped.ix = ir2vec(shift(IRanges(values(out)$query.start, values(out)$query.end), query.offsets), each = x.expand);
                unmapped.ix = setdiff(1:query.length, mapped.ix)
                
                ## look up y vals corresponding to mapped positions
                ## expand these vals when |scale| < 1,
                ## reversing where appropriate (ie scale < 0)                
                df.vals = data.frame(chr = as.character(Rle(as.character(seqnames(out)), width(out))),
                  pos = ir2vec(ranges(out), link.scales<0), stringsAsFactors = F)

                ## if any links scale down
                ## need to expand df.vals by appropriate link.scales, being mindful of left and right edge issues
                down.scaled = abs(link.scales)<1
                if (any(down.scaled)) {
                    ## length should = nrow(df.vals)
                    ## expand represents how many times we plan to copy each row of df.vals
                    y.expand = rep(1, length(link.scales));
                    y.expand[down.scaled] = 1/abs(link.scales)[down.scaled]

                    xq.starts = start(x[values(out)$query.id]) + values(out)$query.start - 1;
                    xq.ends = start(x[values(out)$query.id]) + values(out)$query.end - 1;
                    
                    ## need to take care of "fractional" edge cases on both left and right of each interval
                    ## ie y coordinates for which there are fewer than y.expand[k] x coordinates assigned 
                    right.edge.ix = cumsum(width(out))
                    left.edge.ix = c(0, right.edge.ix[1:(length(right.edge.ix)-1)])[1:length(right.edge.ix)] + 1
                    mod.left = y.expand[down.scaled]-
                      ((xq.starts[down.scaled] - start(.Object@.galx[values(out)$link.id])[down.scaled]) %% y.expand[down.scaled])
                    mod.right = 1+((xq.ends[down.scaled]-start(.Object@.galx[values(out)$link.id])[down.scaled]) %% y.expand[down.scaled])
                    
                    ## expand y.expand to have length = nrow(df.vals)
                    y.expand = as.integer(Rle(y.expand, width(out)));  
                    y.expand[right.edge.ix[down.scaled]] = mod.right
                    y.expand[left.edge.ix[down.scaled]] = mod.lefts
                    
                    ## now replicate the appropriate rows of df.vals
                    df.vals = df.vals[as.integer(Rle(1:length(y.expand), y.expand)), ];
                }

                tmp.out = df.out[mapped.ix, ]
                tmp.out[, c('chr', 'pos')] = df.vals[, c('chr', 'pos')]                
                df.out = rbind(df.out[unmapped.ix, ], tmp.out)
                df.out = df.out[order(df.out$query.coord), ]

                ## df.vals should have the same dimension as ix
                ## use to populate df.out

                if (format == 'df'){
                    df.out = df.out[!duplicated(df.out[, c('query.coord')]), ];
                    if (input.grl){
                        df.out = split(df.out, df.out$list.id)
                    }
                } else {
                    if (input.grl){
                        df.out = split(df.out, df.out$list.id)                    
                    }

                }
                out = df.out;
            } else {                
                if (input.grl){
                    if (length(out)>0){
                            out = split(out, values(out)$list.id)
                            
                            if (!is.null(grl.names)){
                                ix = match(names(out), grl.names)
                                names(out) = grl.names[ix]
                            } else{
                                ix = names(out);                                    
                            }

                            if (ncol(grl.val)>0){
                                values(out) = grl.val[ix, ,drop = FALSE];
                            }
                            
                            if (split.grl){
                                out = grl.split(out)
                            }
                    } else{
                        out = gr.fix(GRangesList(), out)
                    }
                }
            }
            
            ## output gTrack if gTrack was the input
            if (!is.null(x.track)){
                out = gTrack(out);
                formatting(out) = formatting(x.track)
                colormap(out) = colormap(x.track);
            }
            
            return(out)            
})




#' @name dim
#' @title dim
#'
#' dimension of chain is @n x @m
#'
#' @export
setMethod("dim", signature(x = "gChain"), function(x) return(c(x@.n, x@.m)))




#' @name links
#' @title links
#'
#' returns links associated with chain
#'
#' @export
setGeneric('links', function(.Object, ...) standardGeneric('links'))
setMethod("links", signature(.Object = "gChain"), function(.Object){
    return(list(x = .Object@.galx, y = .Object@.galy))                             
})



#' @name values
#' @title values
#' @description
#'
#' Print meta data of gChain
#'
setMethod("values", signature(x= "gChain"), function(x){
    return(x@values)
})




#' @name scale
#' @title scale
#' @description
#'
#' Get scale of gChain
#'
setMethod("scale", signature(x= "gChain"), function(x){
    return(unique(abs(x@.scale)))
})




#' @name pads
#' @title pads
#' @description
#'
#' Get pads of gChain
#'
setGeneric('pads', function(.Object, ...) standardGeneric('pads'))
setMethod("pads", signature(.Object = "gChain"), function(.Object){
    return(data.frame(pad.left = .Object@.pad.left, pad.right = .Object@.pad.right, stringsAsFactors = F))
})




#' @name genomes
#' @title genomes
#' @description
#'
#' return seqinfos associated with domain and range genome of gChain
#'
setGeneric('genomes', function(.Object, ...) standardGeneric('genomes'))
setMethod("genomes", signature(.Object = "gChain"), function(.Object){
    return(lapply(seqinfo(.Object), seqinfo2gr))
})



setMethod('$', 'gChain', function(x, name){
    return(links(x)[[name]])
})




#' @name c
#' @title c
#' @description
#'
#' Concatenate gChains, i.e. concatenate their links
#' they must have the same domain and range genome
#'
#' @export
setMethod('c', 'gChain', function(x, ...){

    if (missing('x')){
        args = list(...)
    } else{
        args = c(list(x), list(...))
    }
    if (any(sapply(args, class) != 'gChain')){
        stop('Error: At least one of the objects to be concatanted is not a gChain')
    }

    return(gChain(unlist(do.call('GRangesList', lapply(args, function(x) x@.galx))),
        unlist(do.call('GRangesList', lapply(args, function(x) x@.galy))),
            pad.left = unlist(lapply(args, function(x) x@.pad.left)),
            pad.right = unlist(lapply(args, function(x) x@.pad.right))))
})



#############################
#' @name expand
#' @title expand
#' @description
#'
#' GChain::expand
#'
#' Expands gChain by "space" base pairs.   By default will pad the x intervals with space and pad and shift the y intervals,
#' however this can be changed with the shift.x and shift.y flags.
#'
#' This is useful when gChain represents a multiple sequence alignment, and we would like to visualize data lifted "around" the
#' alignment, e.g. some window upwards or downwards
#' 
#' If abs(scale) of the gChain != 1, then space will be added to the "smaller" space and space*scale positions will be added to the
#' larger space, i.e. if abs(scale)<1 then y is the smaller space, if abs(scale)>1, then x is the smaller space
#'
#' if space = Inf, then the Chain will map the entire sequence length (ie the whole chromosome)
#' (of x or y, whichever is the largest for the given interval pair)
#'
#############################
setMethod("expand", signature(x = "gChain"), function(x, space = NULL, shift.x = FALSE, shift.y = TRUE){
    
    abs.scale = unique(abs(x@.scale))

    if (is.null(space)){
        space = Inf;
    }
            
    if (abs.scale<=1){
        if (is.infinite(space)){
            slen.x = floor(seqlengths(x@.galx)[as.character(seqnames(x@.galx))]*abs.scale)
            slen.y = seqlengths(x@.galy)[as.character(seqnames(x@.galy))]
            space = pmax(slen.x, slen.y)
        }
        
        right.space = space;
        left.space = pmin(pmin(space, start(x@.galy) + as.numeric(shift.y)*space-1),
                          floor((start(x@.galx) + (1/abs.scale)*as.numeric(shift.x)*space-1)*abs.scale));
        
        tmp = aggregate(end(x@.galy) + right.space + as.numeric(shift.y)*right.space,
                        by = list(as.character(seqnames(x@.galy))), FUN = max)
        new.slen.y = structure(tmp[,2], names = tmp[,1])
        seqlengths(x@.galy) = pmax(seqlengths(x@.galy), new.slen.y[seqlevels(x@.galy)])

        tmp = aggregate(end(x@.galx) + right.space * (1/abs.scale) + as.numeric(shift.x)*right.space*(1/abs.scale),
                        by = list(as.character(seqnames(x@.galx))), FUN = max)
        new.slen.x = structure(tmp[,2], names = tmp[,1])
        seqlengths(x@.galx) = pmax(seqlengths(x@.galx), new.slen.x[seqlevels(x@.galx)])

        if (shift.x){
            x@.galx = GenomicRanges::shift(x@.galx, 1/abs.scale*(right.space))
        }
        
        if (shift.y){
            x@.galy = GenomicRanges::shift(x@.galy, right.space)
        }

        x@.galx = gr.pad(x@.galx, cbind((1/abs.scale)*left.space, (1/abs.scale)*right.space))                
        x@.galy = gr.pad(x@.galy, cbind((left.space), (right.space)))
    } else if (abs.scale>1){
        if (is.infinite(space)){
            slen.x = seqlengths(x@.galx)[as.character(seqnames(x@.galx))]
            slen.y = floor(seqlengths(x@.galy)[as.character(seqnames(x@.galy))]/abs.scale)
            space = pmax(slen.x, slen.y)
        }
              
        right.space = space;
        left.space = pmin(pmin(space, floor((start(x@.galy) + as.numeric(shift.y)*space*abs.scale-1)/abs.scale),
            start(x@.galx) + as.numeric(shift.x)*space-1));
                
        tmp = aggregate(end(x@.galx) + right.space + as.numeric(shift.x)*right.space,
            by = list(as.character(seqnames(x@.galx))), FUN = max)
        new.slen.x = structure(tmp[,2], names = tmp[,1])
        seqlengths(x@.galx) = pmax(seqlengths(x@.galx), new.slen.x[seqlevels(x@.galx)])
                
        tmp = aggregate(end(x@.galy) + right.space*abs.scale + as.numeric(shift.y)*right.space*abs.scale,
            by = list(as.character(seqnames(x@.galy))), FUN = max)
        new.slen.y = structure(tmp[,2], names = tmp[,1])
        seqlengths(x@.galy) = pmax(seqlengths(x@.galy), new.slen.y[seqlevels(x@.galy)])
                

        if (shift.x){
            x@.galx = GenomicRanges::shift(x@.galx, right.space)
        }
    
        if (shift.y){
            x@.galy = GenomicRanges::shift(x@.galy, abs.scale*right.space)
        }

                
        x@.galx = gr.pad(x@.galx, cbind(left.space, right.space))
        x@.galy = gr.pad(x@.galy, cbind(abs.scale*(left.space), abs.scale*(right.space)))
    }              
            
    validObject(x)
    return(x);
})


setReplaceMethod("values", signature(x= "gChain"), function(x, value){
    x@values = value;
    return(x)
})




#' @name seqinfo
#' @title seqinfo
#' @description
#'
#' Return seqinfo for gChain
#'
#' @export
setMethod("seqinfo", signature(x = "gChain"), function(x){
    return(list(x = seqinfo(x@.galx), y = seqinfo(x@.galy)))
})




#' @name t
#' @title t
#' @description
#'
#' "transpose of chain" is flipping x and y variables (i.e. inverting the mapping)
#'
#' @export
setMethod("t", signature(x = "gChain"), function(x){
    tmp.n = x@.n
    tmp.galx = x@.galx

    x@.galx = x@.galy
    x@.n = x@.m

    x@.galy = tmp.galx
    x@.m = tmp.n
    x@.scale = 1/x@.scale;
            
    return(x)
})


#############################
#' @name breaks
#' @title breaks 
#' @description
#'
#' Method to extract "breaks" from a gChain mapping genome x to y
#' i.e. pairs of positions that are contiguous in y but were non-contiguous in x
#'
#' Returns a GRangesList in genome x containing the location of pairs of positions
#' that are contiguous in y that were not contiguous in x.
#' 
#' Breaks are only defined in one direction (ie x to y)
#' (if you want to do reverse then just "transpose" the chain, or use rev = TRUE)
#'
#' @export
#############################
setGeneric('breaks', function(x, ...) standardGeneric('breaks'))
setMethod("breaks", signature(x = "gChain"), function(x, rev = FALSE) {

    if (!inherits(x, 'gChain')){
        stop("Input must be a gChain object.")
    }

    if (rev){
        x = t(x)
    }

    if (length(x@.galx)==0){
        return(GRangesList())
    }
            
    seed = suppressWarnings(c(GenomicRanges::shift(gr.start(x@.galy, width=2, clip=FALSE)-1), GenomicRanges::shift(gr.end(x@.galy, width=2, clip=FALSE), 1)))
    strand(seed) = '+'
    
    ## only lift junctions that have not fallen over the edge (i.e. at the beginning or end of a seq)
    seed = seed[width(seed)>1]
    seed.l = lift(t(x), seed);
              
    ## "broken" ranges will yield at least two ranges when lifted, and will have a width
    ## of 1 grab their ids
    seed.l = seed.l[width(seed.l)==1]            

    if (length(seed.l)==0){
        return(GRangesList())
    }
              
    broken.id = as.numeric(names(which(vaggregate(seed.l$query.start, by = list(seed.l$query.id), FUN = function(x) any(x==1) & any(x==2)))))
    seed.l = seed.l[seed.l$query.id %in% broken.id]
           
    if (length(seed.l)==0){
        return(GRangesList())
    }

    ## find all pairs of lifted ranges that are no longer contigous
    ## ie one range has query.start = 1 and the other query.start = 2
    ##    and they both have the same query.id
    ## and they are on different seqnames in the y genome or more than 1 bp apart
    ##
    qid = seed.l$query.id
    qix = seed.l$query.start
    id = 1:length(seed.l)            
    ij = do.call('rbind', lapply(broken.id, function(x){
        y = which(qid==x); z = qix[y]; zy1 = y[z==1]; zy2 = y[z==2];
        return(cbind(rep(zy1, length(zy2)), rep(zy2, each = length(zy1))))
    }));

    keep = (end(seed.l)[ij[,2]] - start(seed.l)[ij[,1]])!=1 | as.logical(seqnames(seed.l)[ij[,2]] != seqnames(seed.l)[ij[,1]])
    ij = ij[keep, , drop = FALSE]

    ## we flip the first range in the pair so it points backward on x
    ## while the second range points forward (towards the next segment on x)
    tmp.out = c(gr.flipstrand(seed.l[ij[,1]]), seed.l[ij[,2]])
    bk.id = rep(1:nrow(ij), 2)
    out = split(tmp.out, bk.id)

    return(out)
})


#' @name cn
#' @title cn
#' @description
#'
#' A method to extract copy number from a chain mapping genome x to y
#'
#' Returns a GRanges in genome x specifying the number of copies of every interval
#' in x that are contained in y
#'
#' as "breaks" only goes in the "forward" direction
#'
#' @export
setGeneric('cn', function(x, ...) standardGeneric('cn'))
setMethod("cn", signature(x = "gChain"), function(x, rev = FALSE){
    if (rev){
        x = t(x)
    }

    g2 = seqinfo2gr(x@.galy);
    g2.l = lift(t(x), g2)
    out = as(coverage(g2.l), 'GRanges')
    colnames(values(out)) = 'cn'
    return(out)                        

})



#' @name [
#' @title [
#' @description
#' 
#' Subsetting links in a gChain, e.g. using features of the links metadata
#' 
#' @export
setMethod('[', 'gChain', function(x, i){
    x@.galx = x@.galx[i]
    x@.galy = x@.galy[i]
    x@.scale = x@.scale[i]
    x@values = as.data.frame(x@values[i, , drop = F])
    x@.pad.left = x@.pad.left[i]
    x@.pad.right = x@.pad.right[i]

    validObject(x)            
    return(x)
})


#' @name gChain
#' @title gChain
#' @description
#'
#' instantiate a new gChain
#'
#' @export
gChain = function(x = NULL, y = NULL, pad.left = 0, pad.right = 0, scale = NULL, val = data.frame()) new('gChain', x=x, y=y, pad.left = pad.left, pad.right = pad.right, scale = scale, val = val)

######
# Basic gChain synthesizers
#
# Functions that make gChains to represent spliced mappings, pairwise alignment, multiple sequence alignment,
# and cigar strings.
#
#######


#' @name spChain (spliced chain)
#' @title spChain
#' @description
#'
#' Takes (named) GRangesList and outputs gChain mapping to new genome where every
#' GRangesList item is a chromosome (e.g. named after that list item, such as transcript)
#' and the local coordinates are determined the splicingthe ranges in the order (or reverse
#' order as specified) of the ranges.
#'
#' @export
spChain = function(grl, rev = FALSE){

    if (is.null(names(grl))){
        names(grl) = as.character(1:length(grl))
    }

    names(grl) = dedup(names(grl), ' copy ');
    values(grl)$rev = rev;
  

    gr.dt = as.data.table(grl)
    intA = GRanges(as.character(gr.dt$seqnames), IRanges(gr.dt$start, gr.dt$end), gr.dt$strand, seqlengths = seqlengths(grl))

    out.dt = gr.dt[, .(
        chr.A = seqnames,
        start.A = start,
        end.A = end,
        width,
        chr.B = group_name,
        str = ifelse(rep(values(grl)$rev, elementNROWS(grl)), '-', '+')
    )]
  
    out.dt[, start.B := if (length(width)==1) 1 else cumsum(c(1, width[-length(width)])), by = chr.B]
    out.dt[, end.B := start.B + width-1]

    seqlengths = gr.dt[, sum(width), by = group_name][, structure(V1, names = group_name)]
    intB = GRanges(out.dt$chr.B, IRanges(out.dt$start.B, out.dt$end.B), seqlengths = seqlengths, strand = out.dt$str)
  
    return(gChain(intA, intB))
} 




#' @name paChain
#' @title paChain
#' @description
#' 
#' Make chain from pairwiseAlignment object
#' 
#' Each BSStringSet must have named sequences (otherwise names will be computed
#' from (perfectly unique) sequences.  Each name can only map to a single sequence
#' 
#' NOTE: if XStringSet objects are not provided, then the input is assumed to be DNA
#' (if assume.dna = T), otherwise AA (if assume.dna = F)
#'
#' Pre-computed pa object can be also provided, but then the sequences from which the
#' pa object was computed should be given as well.
#'
#' Returns gChain with seqinfos corresponding to unique seqnames in seq1 and seq2
#' and the interval pairs mapping those two sequences
#' 
#' @export
#' @author Marcin Imielinski
paChain = function(seq1, seq2,
    sn1 = NULL, ## character vector of seqnames for seq1, should  be same length as seq1
    sn2 = NULL, ## character vector of seqnames for seq2, should  be same length as seq2
    gr1 = NULL, ## GRanges corresponding to seq1, MUST be of same width, length, and names as seq1
    gr2 = NULL, ## GRanges corresponding to seq2, MUST be of same width, length, and names as seq2  
    pa = NULL, 
    verbose = TRUE,
    sl1 = NULL, ## seqlengths vector for seq1 (useful if provided seqs exist in larger set)
    sl2 = NULL, ## seqlengths vector for seq2 (useful if provdied seqs if exist in larger set)  
    score.thresh = NULL, ## optional numeric scalar specifying lowest score alignment to include in chain
    both.strands = FALSE, # both strands = T will try to align seq1 to both strands of seq2
    keep.best = TRUE, # if both.strands = T, keep.best = T will only keep the best seq2 alignment (forward vs backward) for each seq1
    assume.dna = TRUE, 
    mc.cores = 1, ## how many cores
    mc.chunks = mc.cores,  ## how many chunks to split the data when farming out to each core
    pintersect = FALSE, ## calls to gr.findoverlaps (and to lift) use pintersect. See gr.findoverlaps. Best for small genomes (e.g. reads)
    ...)
{

    if (inherits(seq1, 'factor')){
        seq1 = as.character(seq1)
    }

    if (inherits(seq2, 'factor')){
        seq2 = as.character(seq2)
    }

    if (inherits(seq1, 'character')){
        if (any(ix <- is.na(seq1))){
            seq1[ix] = 'NA'
        }
                
        if (assume.dna){
            seq1 = DNAStringSet(seq1)
        } else{
            seq1 = AAStringSet(seq1)
        }
    }
  
    if (inherits(seq2, 'character')){
        if (any(ix <- is.na(seq2))){
            seq2[ix] = 'NA'
        }
      
        if (assume.dna){
            seq2 = DNAStringSet(seq2)
        } else{
            seq2 = AAStringSet(seq2)
        }
    }

    if (inherits(seq1, 'DNAString')){
        seq1 = DNAStringSet(seq1)
    }
    if (inherits(seq1, 'AAString')){
        seq1 = AAStringSet(seq1)
    }
  
    if (inherits(seq2, 'DNAString')){
        seq2 = DNAStringSet(seq2)
    }
    if (inherits(seq2, 'AAString')){
        seq2 = AAStringSet(seq2)
    }
  
    if (!(inherits(seq1, 'XStringSet')) | !(inherits(seq2, 'XStringSet'))){
        stop('seq1 and seq2 arguments must be both XStringSet objects (e.g. DNAStringSet, RNAStringSet, AAStringSet).')
    }

    if (is.null(sn1)){
        if (is.null(names(seq1))){
            sn1 = as.character(1:length(seq1))
        } else{
            sn1 = names(seq1)
        }
    } else{
        sn1 = cbind(1:length(seq1), as.character(sn1))[,2]
    }
  
    if (is.null(sn2)){
        if (is.null(names(seq2))){
            sn2 = as.character(1:length(seq2))
        }else{
            sn2 = names(seq2)
        }
    } else{
        sn2 = cbind(1:length(seq2), as.character(sn2))[,2]
    }
   
    ## if sequence lengths are not equal, then try to replicate shorter
    if (length(seq1) != length(seq2)){
        length.rat = length(seq1)/length(seq2)
      
        if ((length.rat %% 1)!=0 & ((1/length.rat) %% 1)!=0){
            stop('seq1 and seq2 must be either of same lengths or have lengths such that the length of the longer object as integer multiple of the length of the shorter, so that the shorter object may be replicated')
        }
      
        if (length.rat>1){
            seq2 = rep(seq2, length.rat)
            sn2 = rep(sn2, length.rat)

            if (!is.null(gr2)){
                gr2 = rep(gr2, length.rat)
            }
        } else{
            seq1 = rep(seq1, 1/length.rat)
            sn1 = rep(sn1, 1/length.rat)

            if (!is.null(gr1)){
                gr1 = rep(gr1, 1/length.rat)
            }
        }      
    }

    if ((width(seq1)==0) | (width(seq2)==0)){
        stop("Either seq1 and/or seq2 is of width 0. This will result in an error with pairwiseAlignment(seq1, seq2).")
    }
  
    ## check to make sure that each sequence name maps to exactly one sequence in both seq1 and seq2    
    if (any(duplicated(sn1)) && FALSE){
        useq1 = unique(as.character(seq1))
        seq1.ix = match(as.character(seq1), useq1)
        uname1 = unique(sn1)
    }
  
    if (any(duplicated(sn2)) && FALSE ){      
        useq2 = unique(as.character(seq2))
        seq2.ix = match(as.character(seq2), useq2)      
        uname2 = unique(sn2)

    }
  
    if (!is.null(gr1)){
        if (length(gr1) != length(seq1)){
            stop('Length of gr1 and seq1 not compatible ')
        }

        if (any(width(gr1) != width(seq1))){
            stop('Widths of gr1 and seq1 not compatible ')
        }

        if (!is.null(names(gr1))){
            names(gr1) = sn1
        }

        if (any(names(gr1) != sn1)){
            stop('Names of gr1 and seq not compatible')
        }    
    }
 
    if (!is.null(gr2)){
        if (length(gr2) != length(seq1)){
            stop('Length of gr2 and seq1 not compatible ')
        }

        if (any(width(gr2) != width(seq1))){
            stop('Widths of gr2 and seq1 not compatible ')
        }

        if (!is.null(names(gr2))){
            names(gr2) = sn2
        }
      
        if (any(names(gr2) != sn2)){
            stop('Names of gr2 and seq not compatible')
        }
    }
      
    if (verbose) {
        time <- proc.time()
        message('Starting alignment\n')
    }
  
  
    if (is.null(pa)){
        if (both.strands){
            uname2 = unique(sn2)
            rc_map = structure(rep(uname2, 2), names = c(uname2, paste(uname2, '_rc', sep = '')))
            rc_strand = structure(rep(c('+', '-'), each = length(uname2)), names = names(rc_map))
            seq1 = rep(seq1, 2)
            sn1 = rep(sn1, 2)
            tmp = structure(reverseComplement(seq2))
            sn2.tmp = paste(sn2, '_rc', sep = '')
            seq2 = c(seq2, tmp);
            sn2 = c(sn2, sn2.tmp)
        }

        opts <- list() ### list(...)
        if ('numchunk' %in% names(opts)){
            numchunk <- opts$numchunk
        } else{
            numchunk <- NA
        }

        if (mc.cores > 1){
            pa <- mc.pairwiseAlignment(seq1, seq2, numchunk=numchunk, mc.cores=mc.cores) ## , ...)
        } else{
            pa <- pairwiseAlignment(seq1, seq2) 
        }
    }

  
    if (verbose) {
        print(proc.time() - time)
        message('Finished aligning\nPrepping chain ..\n')
        time <- proc.time()
    }

    uname1 = unique(sn1);
    uname2 = unique(sn2);

    uix.1 = match(uname1, sn1)
    uix.2 = match(uname2, sn2);
  
    ## Now build two gChains mapping the two sequence spaces to the "alignment space"
    ## where each alignment yields its own chromosome
    ## define alignment, subject, and pattern "genomes"

    ## if both.strands and keep.best we choose only a single alignment for each rev comp pair
    if (both.strands & keep.best){
        sc = score(pa)

        ## new way with data.table. 21 seconds on length(sn1) ~ 6 million
        dt <- data.table(id=seq_along(seq2), sn1=sn1, rc.map=rc_map[sn2], sc=sc, key=c('sn1', 'rc.map'))
        setkey(dt, sn1, rc.map, sc) # sort by sn1, then rc.map, then sc

        ## hack because fromLast not supported for data.table
        isdup <- duplicated(paste(dt$sn1, dt$rc.map), fromLast=TRUE)
        #isdup <- duplicated(dt, by=c('sn1', 'rc.map'))
        keep.ix <- dt$id[!isdup]
        ## old way. 380 seconds on same data. After sorting dropping names, identical=TRUE
        #system.time(keep.ix.df <- vaggregate(1:length(seq2), by = list(sn1, rc_map[sn2]), function(x) x[which.max(sc[x])]))
    } else{
        keep.ix = 1:length(seq1);
    }

    sct <- score(pa) > score.thresh
    if (!is.null(score.thresh)){
        keep.ix = intersect(keep.ix, which(sct))
    }

    pa = pa[keep.ix];        
   
    if (is.null(sl1)) {
        system.time(sl1 <- width(seq1)[uix.1])
        names(sl1) <- sn1[uix.1]
    } else {
        ## warning. intractably slow for large uix.1, etc
        sl1[sn1[uix.1]] = width(seq1)[uix.1];   
    }

    if (is.null(sl2)) {
        system.time(sl2 <- width(seq2)[uix.2])
        names(sl2) <- sn2[uix.2]
    } else {
        ## warning. intractably slow for large uix.2, etc
        sl2[sn2[uix.2]] = width(seq2)[uix.2];   
    }

    if (verbose){
        print(proc.time() - time)
    }

    if (length(pa)==0){
        return(gChain(GRanges(seqlengths = sl1), GRanges(seqlengths = sl2)))
    }

    sinfo1 = Seqinfo(seqlengths = sl1, seqnames = names(sl1))
    sinfo2 = Seqinfo(seqlengths = sl2, seqnames = names(sl2))
    sinfo.ali = Seqinfo(seqlengths = nchar(as.character(pa)), seqnames = as.character(1:length(pa)))

    ## deletions specify pattern gaps in alignment coordinates
    del.len = elementNROWS(deletion(pa))
    del.ix = unlist(sapply(1:length(del.len), function(x) rep(x, del.len[x])))

    if (length(del.ix)>0){
        ## each seqname represents a different index of pa
        gr.del = GRanges(as.character(del.ix), unlist(deletion(pa)), seqlengths = seqlengths(sinfo.ali), strand = '+')
        ## FIX to Biostrings pairwise alignment bug / issue
        gr.del = GenomicRanges::shift(gr.del, levapply(width(gr.del), as.character(seqnames(gr.del)), function(x) if (length(x)>1) cumsum(c(0, x[1:(length(x)-1)])) else return(0)))
    } else{
        gr.del = GRanges(seqlengths = seqlengths(sinfo.ali)) 
    }

    pat.mapped = setdiff(seqinfo2gr(sinfo.ali), gr.del) ## pattern corresponds to all non-deleted positions in the alignment    
    tmp = spChain(split(pat.mapped, seqnames(pat.mapped))) ## map to these to pattern coordinates
    al.ix = as.numeric(as.character(seqnames(links(tmp)$y)))
    ali2pat = gChain(links(tmp)$x, GRanges(sn1[keep.ix[al.ix]], ## make sure to shift to account for start gap
        ranges = GenomicRanges::shift(ranges(links(tmp)$y), start(Biostrings::pattern(pa[al.ix]))-1), strand = '+', seqlengths = seqlengths(sinfo1)))
  
    ## insertions specify subject gaps in alignment coordinates
    ins.len = elementNROWS(insertion(pa))
    ins.ix = unlist(sapply(1:length(ins.len), function(x) rep(x, ins.len[x])))
  
    if (length(ins.ix)>0){
        ## each seqname represents a different index of pa
        gr.ins = GRanges(as.character(ins.ix), unlist(insertion(pa)), seqlengths = seqlengths(sinfo.ali), strand = '+')
        ## FIX to Biostrings pairwise alignment bug
        gr.ins = GenomicRanges::shift(gr.ins, levapply(width(gr.ins), as.character(seqnames(gr.ins)), function(x) if (length(x)>1) cumsum(c(0, x[1:(length(x)-1)])) else return(0)))
    } else{
        gr.ins = GRanges(seqlengths = seqlengths(sinfo.ali))
    }
  
    subj.mapped = setdiff(seqinfo2gr(sinfo.ali), gr.ins) ## subject corresponds to all non-inserted positions in alignment coordinates
    tmp = spChain(split(subj.mapped, seqnames(subj.mapped))) ## map these to subject coordinates
    al.ix = as.numeric(as.character((seqnames(links(tmp)$y)))) ### have to do as.numeric(as.vector( here .. because of weird factor issues
    ali2subj = gChain(links(tmp)$x, GRanges(sn2[keep.ix[al.ix]],
        ranges = GenomicRanges::shift(ranges(links(tmp)$y), start(subject(pa[al.ix]))-1),  ## make sure to shift to account for start gap
        strand = '+', seqlengths = seqlengths(sinfo2)), val = data.frame(score = score(pa[as.numeric(as.character(seqnames(subj.mapped)))])))
  
    if (verbose) {
        print(proc.time() - time)
        message('Prepped chain\nMultiplying\n')
        time <- proc.time()
    }

    # this new chain maps pat-alignment combinations to subj-alignment combinations
    # now need to merge seqnames on x and y links so that we are mapping
    # subject to pattern
  
    out <- ali2subj * t(ali2pat)

    if (verbose) {
        print(proc.time() - time)
        message('Finalizing chain\n')
        time <- proc.time()
    }
    if (verbose & both.strands){
        message('Strand collapsing\n')
    }

    if (both.strands){
        sl2 = seqlengths(sinfo2)
        tmp = unique(rc_map[names(sl2)]);
        sl2.2 = structure(sl2[tmp], names = tmp)
        tmp.gr1 = GRanges(names(sl2), IRanges(1, sl2), strand = '+', seqlengths = sl2)      
        tmp.gr2 = GRanges(rc_map[names(sl2)], IRanges(1, sl2), strand = rc_strand[names(sl2)], seqlengths = sl2.2)
        out = gChain(tmp.gr1, tmp.gr2) * out
    }

    ## lift back through gr1 
    if (!is.null(gr1)){
        if (verbose){
            message('Lifting through gr1\n')
        }
        si1 = seqinfo2gr(links(out)$x)
        out = out * gChain(gr1[seqnames(si1)], si1)
    }

    ## lift back through gr1 
    if (!is.null(gr2)){
        if (verbose){
            message('Lifting through gr1\n')
        }
        si2 = seqinfo2gr(links(out)$x)
        out = gChain(si2, gr1[seqnames(si1)]) * out
    }

    if (verbose){
        print(proc.time() - time)
    }

  return(out)

}


#' @name cgChain
#' @title cgChain
#' @description
#' cgChain (i.e. "CIGAR")
#'
#' Processes a pairwise alignment from read to genomic coordinates inputted as a vector of CIGAR strings
#' representing edit operations from subject to pattern.
#' ie an insertion corresponds a gap in the subject and deletion corresponds to a gap in the pattern
#' returns a chain that maps subject to pattern (e.g. genome space to read space)
#'
#' Also, "cigar" can be a GRanges or GappedAlignment object with $cigar field where the ranges specify flanking subject coordinates
#' of the aligned pattern (e.g. read).  As by convention, subject coordinates are assumed to NOT contain flanking soft clips.
#' Internal soft clips (if any exist - this would be unconventional use) are treated as insertions.
#' If seqnames ("sn") input not provided then the read seqnames will be the $qname field + IsFirstMateRead part of read flag.
#'
#' Returns a gChain mapping read space to genome space.
#' cigar --- ## can be character vector or GRanges GAlignment with $cigar and $qname string field)
#' 
#' @export
#' @author Marcin Imielinski
cgChain = function(cigar, sn = NULL, verbose = TRUE){
    
    gr = NULL;
    cig.names = NULL;

    if (inherits(cigar, 'GappedAlignment')){
        if (any(ix <- is.na(values(cigar)$cigar))){
            warning('Some CIGAR values are NA, ignoring these ranges.')
            cigar = cigar[!ix]
        }
        gr = granges(cigar)

        if (!is.null(sn)){
            cig.names = cbind(1:length(cigar), sn)[,2]
        } else{
            cig.names = paste(values(cigar)$qname, ifelse(bamUtils::bamflag(values(cigar)$flag)[, 'isFirstMateRead'], 1, 2))
        } 
        
        cigar = structure(values(cigar)$cigar, names = cig.names)

    } else if (inherits(cigar, 'GRanges')) {
        if (any(ix <- is.na(values(cigar)$cigar))){
            warning('Some CIGAR values are NA, ignoring these ranges.')
            cigar = cigar[!ix]
        }
        gr = cigar;

        if (!is.null(sn)){
            cig.names = cbind(1:length(cigar), sn)[,2]
        } else{
            cig.names = paste(values(cigar)$qname, ifelse(bamUtils::bamflag(values(cigar)$flag)[, 'isFirstMateRead'],'1',
                ifelse(bamUtils::bamflag(values(cigar)$flag)[, 'isSecondMateRead'], '2', '')), sep = '')
        }
        
        cigar = structure(cigar$cigar, names = cig.names)
    } else if (inherits(cigar, 'data.table')){
        ix <- !is.na(cigar$cigar)
        cigar <- cigar[ix]
        gr <- cigar
        gr = dt2gr(gr)
        if (!is.null(sn)){
            cig.names <- sn[ix] #cbind(seq_along(nrow(cigar)), sn)[,2]
        } else{
            cig.names <- cigar$qname
        }
        cigar <- structure(cigar$cigar, names=cig.names)
    } else{
        if (!is.null(sn)){
            cig.names = cbind(1:length(cigar), sn)[,2]
        } else{
            cig.names = 1:length(cigar)
        } 
    }
              
    if (is.null(cig.names) || "" %in% cig.names){
        cig.names = as.character(1:length(cigar))
    }

    if (verbose){
        message('Parsing CIGARs\n')
    }

    #suppressWarnings(cigar.s <- explodesplitCigar(cigar))
    #suppressWarnings(cigar.s <- CigarOpLengths(cigar))
    cigar.m <- GenomicAlignments::explodeCigarOps(cigar)
    cigar.l <- GenomicAlignments::explodeCigarOpLengths(cigar)
    
    ## reverse CIGARs for ranges mapped to the negative strand
    if (!is.null(gr))
        if (any(ix <- as.logical(strand(gr)=='-'))){
            cigar.l[ix] = lapply(cigar.l[ix], rev)
            cigar.m[ix] = lapply(cigar.m[ix], rev) 
    }

    ## D = charToRaw('D')
    ## S = charToRaw('S')
    ## H = charToRaw('H')
    ## I = charToRaw('I')

    D = 'D'
    S = 'S'
    H = 'H'
    I = 'I'
    
    sinfo.cigar <- Seqinfo(seqlengths = sapply(cigar.l, sum), seqnames = as.character(1:length(cigar.l))) ## 3.1
    isnotd = lapply(cigar.m, function(x) x != D) ## 3.1
    sl.pat <- vaggregate(sapply(seq(length(isnotd)), function(x) sum(cigar.l[[x]][isnotd[[x]]])), by = list(cig.names), FUN = max, na.rm = T) ## 3.1
       
    isnotsi = lapply(cigar.m, function(x) !(x %in% c(I,S,H))) ## 3.1
    sinfo.ref <- Seqinfo(seqlengths = sapply(seq(length(isnotsi)), function(x) sum(cigar.l[[x]][isnotsi[[x]]])), seqnames = as.character(1:length(cigar))) # 3.1
    sinfo.template <- Seqinfo(seqlengths = sl.pat, seqnames = names(sl.pat))        

    ## digest / unlist cigar some more
    cig.id    <- unlist(lapply(seq_along(cigar.l), function(x) rep(x, length(cigar.l[[x]]))))
    cig.type  <- unlist(cigar.m)
    cig.wid   <- unlist(cigar.l) 
    cig.start <- levapply(cig.wid, cig.id, function(x) if (length(x)>1) cumsum(c(1, x[1:(length(x)-1)])) else return(1))

    if (verbose){
        message('Preparing gChains\n')    
    }

    ##
    ## now we need to map template --> CiGARs --> reference coordinates
    ## by building chains and multiplying them 
    ##
    
    ## prep ranges corresponding to insertion and deletion gaps in alignment space
    ## by only taking into account gaps we should  create a cleaner chain, i.e. with fewer
    ## segments

    ## D are gaps in the template        
    ix = cig.type == D
    ## populate with template gaps
    if (any(ix)){
        gr.del = GRanges(cig.id[ix], IRanges(cig.start[ix], width = cig.wid[ix]), seqlengths = seqlengths(sinfo.cigar), strand = '+')
    } else{
        ## otherwise empty granges
        gr.del = GRanges(seqlengths = seqlengths(sinfo.cigar))
    }

    ## ISH are gaps in the reference
    ix = cig.type %in% c(I,S,H)
    if (any(ix)){
        ## populate with reference gaps
        gr.ins = GRanges(cig.id[ix], IRanges(cig.start[ix], width = cig.wid[ix]), seqlengths = seqlengths(sinfo.cigar), strand = '+')
    } else{
        ## empty granges
        gr.ins = GRanges(seqlengths = seqlengths(sinfo.cigar))
    }
        
    sigr.cigar <- seqinfo2gr(sinfo.cigar) ## coordinate space of the CIGAR
    
    ## cigar2template maps cigar to template
    template.mapped <- setdiff(sigr.cigar, gr.del) ## template maps to all non-deleted positions in the CIGAR
    sp  <- split(template.mapped, seqnames(template.mapped))
    tmp <- spChain(sp)
    cigar2template <- gChain(links(tmp)$x,  GRanges(cig.names[as.numeric(as.character(seqnames(links(tmp)$y)))], ranges(links(tmp)$y), strand = '+', seqlengths = seqlengths(sinfo.template)))

    ## cigar2ref maps cigar to reference
    ref.mapped <- setdiff(sigr.cigar, gr.ins) ## ref maps to all non-inserted and non-clipped positions in the CIGAR
    tmp = spChain(split(ref.mapped, seqnames(ref.mapped)))
    cigar2ref = gChain(links(tmp)$x, GRanges(seqnames(links(tmp)$y),
      ranges(links(tmp)$y), strand = '+', seqlengths = seqlengths(sinfo.ref)))
    
    if (verbose){
        message('Multiplying gChains\n')
    }
    
    out <- cigar2template*t(cigar2ref)

    if (verbose){
        message('Finalizing gChain\n')
    }

    ## if CIGAR provided as GRange on "real genome" then lift sinfo.ref onto the genome
    ## keepin gtrack of strand 
    if (!is.null(gr)) {
        ## harmonize gr with cigar .. sometimes right side not consistent
        width(gr) = sapply(1:length(cigar.l), function(x) sum(cigar.l[[x]][cigar.m[[x]] %in% c("D", "M")]))        
        out = out * gChain(gr, GRanges(seqnames(sinfo.ref), IRanges(1, width = width(gr)), strand = '+'), val = as.data.frame(values(gr)))
    }
    
    return(t(out))
}



#' @name maChain
#' @title maChain
#' @description
#'
#' Multiple sequence alignment chain
#'
#' Takes input GRangesList of interval coordinates corresponding to sequences
#' contributing to alignments in pali. 'pali' is a list of alignment matrices
#' (containing letters from the original sequences and '-', '.' for gaps)
#'
#' 'grl' and 'pali', have the same names (corresponding to the alignment, e.g. a protein domain)
#' Also, this property holds:
#' sapply(pali, nrow) = sapply(grl, length)
#'
#' trim will remove all alignment columns that have only gaps
#'
#' alternatively, grl can be an XStringSet object with a "names" attribute containing
#' Biostrings (eg AA or DNA) and gap ("-") characters.  Can also be a list of XStringSet objects.
#'
#' if grl is NULL, then it is assumed that the first (last) non gapped character in each MSA alignment entry
#' represents the first (last) item of the sequence 
#'
#' @export
#' @author Marcin Imielinski
maChain = function(grl = NULL, pali, pad = 0, trim = TRUE, trim.thresh = 0 ### number between 0 and 1 specifying what percentage of alignments a position need to be present in in order to be retained in output coordinates
  )
{
    ## will only make chain for a single alignment
    if (inherits(pali, 'XStringSet') | inherits(pali, 'XString')){
        pali = list(pali)
    } 

    if (is.null(names(pali))){
        names(pali) = 1:length(pali)
    }
    ## will
    if (inherits(pali[[1]], 'XStringSet') | inherits(pali[[1]], 'XString') ){
        npali = names(pali)
        
        pali = lapply(pali, function(this.pali) t(matrix(unlist(strsplit(as.character(this.pali), '')), ncol = length(this.pali), dimnames = list(NULL, names(this.pali)))))
        names(pali) = npali;
    }
    ## assume that pali rows represent the entire sequence
    if (is.null(grl)){
        ## Error in validObject(.Object) : 
        ##   invalid class “GRanges” object: 'mcols(x)' is not parallel to 'x'

        grl = do.call('GRangesList', lapply(pali, function(this.pali) GRanges(rownames(as.data.frame(this.pali)), IRanges(start = 1, rowSums(as.data.frame(this.pali) != '-')))))
        names(grl) = names(pali)
    }
           
    if (!inherits(grl, 'GRangesList')){
        grl = GRangesList(grl)
    }

    if (!inherits(pali, 'list')){
        pali = list(pali)
        names(pali) = 1:length(pali)
    }

    if (is.null(names(grl))){
        names(grl) = names(pali)
    }
    
    
    if (!identical(names(grl), names(pali))){
        stop('Names of grl and pali should be identical')
    }

    if (trim){
        pali = lapply(pali, function(x){
            ix = colSums(x!='.' & x!='-')/nrow(x)>trim.thresh
            x[, ix, drop = FALSE]
        })
    }
    
    slen = lapply(pali, ncol)
    pali.irl = lapply(pali, function(y) apply(y, 1, function(x) as(x != '-' & x != '.', 'IRanges')))
    grl.irl = unlist(lapply(pali.irl, function(y) lapply(y, squeeze)))
    grl.ir = do.call('c', grl.irl)
    pali.ir = unlist(pali.irl)
    gr.ix = unlist(lapply(1:length(pali.ir), function(x) rep(x, length(pali.ir[[x]]))))
    pali.ix = unlist(lapply(1:length(pali.irl), function(x) rep(x, sum(sapply(pali.irl[[x]], length)))))

    gr = unlist(grl);

    intA = GRanges(seqnames(gr)[gr.ix], IRanges::shift(unlist(IRangesList(grl.ir)), start(gr)[gr.ix]-1),
                   strand = strand(gr)[gr.ix], seqlengths = seqlengths(gr));

    intB = GRanges(names(pali)[pali.ix], unlist(IRangesList(pali.ir)), strand = '+', seqlengths = sapply(pali, ncol))

    if (pad==0){
        return(gChain(intA, intB))
    } else{
        return(expand(gChain(intA, intB), pad))        
    }
}


#' @name txChain
#' @title txChain
#' @description
#' 
#' Transcript chain
#'
#' Takes a GRangesList on genomic coordinates (i.e. exons comprising transcripts)
#' and creates a gChain representing mapping onto transcript (or translated protein) coordinates
#'
#' @author Marcin Imielinski
#' @export
txChain = function(grl, txname = NULL, translate = FALSE, exonFrame.field = 'exon_frame',
  val = NULL ## data frame of nrow = length(grl) specifying val to add to chain
  )
{
    if (is.null(txname)){
        txname = names(grl)
    }
    
    if (!is.null(val)){
        if (!is.data.frame(val)){
            val = as.data.frame(val)
        }

        if (is.null(nrow(val)) | (nrow(val) != length(grl))){
            stop('Error: val data.frame must correspond to grl, i.e. have nrows = length(grl)')
        }
    }

    if (!is.null(txname)){
        if (any(ix <- is.na(txname))){
            warning('Ignoring ranges with NA txnames')
            grl = grl[!ix]
            txname = txname[!ix]
          
            if (!is.null(val)){
                val = val[!ix, , drop = FALSE]               
            }
        }
    } else{
        txname = as.character(1:length(grl))
    }
            
    names(grl) = NULL;
    gr = grl.unlist(grl);
    ix = order(c(-1, 1)[1+as.numeric(strand(gr)=='+')] * start(gr))
    gr = gr[ix]
    
    ## in case of refgene we take into account "exon frame" field to correct the frame of exons
    if (exonFrame.field %in% names(values(gr))){
        .fix_frame = function(x) {
            if (length(x)>1){
                g = rep(0, length(x));
                p = rep(1, length(x));
                for (i in 2:length(x)){
                    p[i] = p[i-1] + w[x[i-1]]
                    g[i-1] = (ef[x[i]]-((p[i]-1) %% 3)) %% 3
                    p[i] = p[i] + g[i-1]
                }
                return(g)
            } else{
                return(0)
            }
        }
        
        w = width(gr)
        ef = values(gr)[, exonFrame.field]
        g = levapply(1:length(gr), gr$grl.ix, .fix_frame)

    } else{
        g = 0
    }


    if (!is.null(val)){
        val = val[gr$grl.ix, , drop = FALSE]
    }

    ## Marcin: fixed issue where refactoring . mapping would behave weirdly if some txnames were duplicated
    ## now, the mapping is purely a function of the grl structure, but then seqnames are assigned according
    ## to txname
    gr2.tmp = gr.stripstrand(gr.refactor(gr, gr$grl.ix, gap = g))
    sl2.tmp = seqlengths(gr2.tmp)
    sl2 = vaggregate(sl2.tmp, by = list(txname[as.numeric(names(sl2.tmp))]), FUN = max)    
    gr2 = GRanges(txname[as.numeric(as.character(seqnames(gr2.tmp)))], ranges(gr2.tmp), strand = strand(gr2.tmp), seqlengths = sl2)
    
    gc = gChain(gr, gr2, val = val)

    if (translate){
        si.gr = seqinfo2gr(seqinfo(gc)[[2]])

        if (any((width(si.gr) %% 3) != 0)){
            warning('Widths of some transcripts are not a multiple of 3.  Truncating widths to mod 3.  Make sure that inputs are valid CDSs.')
            end(si.gr) = floor(end(si.gr)/3)*3
        }
        
        gc2 = gChain(si.gr, GRanges(seqnames(si.gr), IRanges(1, end(si.gr)/3), strand = '*', seqlengths = structure(width(si.gr)/3, names = as.character(seqnames(si.gr)))))
        gc = gc2*gc;        
    }

    return(gc)

}


#######
# Rearrangement chain synthesizers
#
# functions to make chains that copy and delete
#
#######

############################################
# duplicate
#
# tandem duplicates a set of gRanges given provided multiplicity
#
# inputs must be nonoverlapping
############################################
duplicate = function(gr, mult = 1, dup.string = 'copy'){
    if (sum(as.numeric(width(reduce(gr)))) != sum(as.numeric(width(gr)))){
        stop('Error: Input ranges must be nonoverlapping')
    }
    return(copy(gr, gr.start(gr), mult = mult, dup.string = dup.string))
}




############################################
# copy
#
# given "from" and "to" GRanges defined on genome A
# returns a gChain from genome A to a second genome B
# which has the same chromosomes but with augmented lengths
# corresponding to copied intervals.
#
# if "to" is a GRanges then the corresponding interval in "from" will be mapped
# to the beginning of the interval in "to".  If "to" is of nonzero width, then
# the interval inside to will be unmapped in the new genome (ie implicitly deleted)
#
# if "to" is left null, then neo-chromosomes are created containing the
# the copied intervals, with names that are derived from the corresponding
# "from" ranges (using gr.tostring)
#
# "to" can be a character vector specifying neochromosomes to copy intervals into ..
# this must be either length 1 or length(from).  If "to" is length 1 or has repeated
# items, then the "from" intervals will be copied (in tandem) to the same target neo chromosome
# If "to" is a character and contains already existing seqlevels, then a warning will be thrown
# up and the "to" chromosome will be renamed with a suffix to make sure it does not conflict
#
# 
############################################
copy = function(from, ## granges of source intervals
        to = NULL, ## granges of target intervals to copy into, or characters vector specifying 'neochromosomes' to copy into ## if null, then names of neochromosomes will be automatically created
        mult = 1, dup.string = 'copy')
{
    genomeA = seqinfo(from)
    old.lens = structure(seqlengths(genomeA), names = seqnames(genomeA))

    if (is.null(to)){
        to = dedup(gr.tostring(from), dup.string)
    }

    if (any(strand(from) == '*')){
        warning('converting some * strands in "from" to +')
        strand(from)[which(as.logical(strand(from)=='*'))] = '+';
    }

    if (is.character(to)){
        values(from)$chr.name = to;
        values(from)$mult = mult;

        # replicate args with respect to multiplicity
        rep.ix = as.integer(Rle(1:length(from), abs(values(from)$mult)))
        from = from[rep.ix]
        to = to[rep.ix]
        values(from)$mult = sign(values(from)$mult)

        # compute lengths of new chroms
        tmp = aggregate(width(from), by = list(values(from)$chr.name), sum)
        tmp = tmp[match(unique(to), tmp[,1]),]

        uchr = unique(values(from)$chr.name)
        chr.ix = lapply(uchr, function(x) which(values(from)$chr.name == x))
        starts.l = lapply(chr.ix, function(x){
            if (length(x)==1){
                1
            } else{
                c(1, 1+cumsum(width(from)[x[2:length(x)]]))
            }
        })
        ends.l = lapply(chr.ix, function(x) cumsum(width(from)[x] + c(0, rep(1, length(x)-1))))
        starts = rep(NA, length(from))
        ends = rep(NA, length(from))
        starts[unlist(chr.ix)] = unlist(starts.l)
        ends[unlist(chr.ix)] = unlist(ends.l)

        new.genome = structure(c(seqlengths(genomeA), tmp[,2]), names = dedup(c(seqnames(genomeA), tmp[,1]), dup.string))
        genomeB = Seqinfo(seqnames = names(new.genome), seqlengths = new.genome);
        
        intA = c(GRanges(seqnames(genomeA), 
            IRanges(rep(1, length(seqlengths(genomeA))), seqlengths(genomeA)), seqlengths = seqlengths(genomeA), strand = '+'),
            GRanges(seqnames(from), ranges(from), seqlengths = seqlengths(genomeA), strand = strand(from)))
        
        intB = c(GRanges(seqnames(genomeA),
            IRanges(rep(1, length(seqlengths(genomeA))), seqlengths(genomeA)), seqlengths = seqlengths(genomeB), strand = '+'),
            GRanges(values(from)$chr.name, IRanges(starts, ends),
                  seqlengths = seqlengths(genomeB), strand = c('-', '+')[1+as.numeric(sign(values(from)$mult>0))]));

        return(gChain(intA, intB, val = data.frame(flag = c(rep(FALSE, length(genomeA)), rep(TRUE, length(from))))))
    } else {
        if (any(strand(to) == '*')){
            warning('converting some * strands in "to" to +')
            strand(to)[which(strand(to)=='*')] = '+';
        }
           
        ## vectorize
        if (length(from) == 1){
            from = rep(from, length(to))
        } else if (length(to) == 1){
            to = rep(to, length(from))
        }
          
        values(from)$mult = mult;
        
        ## find overlapping "to" locations and stack any intervals that map to a given "to"
        ## locus .. note this stacking takes care of an ambiguity arising if the user decides
        ## to map many "froms" to the same or overlapping "to"
        ## ie any nonzero width of each (reduced) "to" interval will be deleted from genomeB
        ## and all "from" intervals will be stacked beginning at the (image) of the start of that
        ## "to" interval on the host genome
        ## stacking will additionally need be done in a chromosome to chromosome fashion ..
        ##        
        ## there are two issues
        ## (1) mapping to the same chromosome, and (2) "to" locations that intersect
        ##
        ## for (1) just have to keep track of widths of all ranges mapping to a single chromosome
        ## and based on their order in the argument and order on that chromosome with respect to start
        ## position you assign them a new coordinate on the new chromosome
        ##
        ## (2) is actually taken care of by (1) except for the situation when we have the start location
        ## for one event located inside the "to" interval for another event .. one argument here would be
        ## that we just throw that event out (since its target site has been deleted by another event)
        ## The other would be that we let those overlapping events stack as if they had the same "start"
        ## and "end" position ... these both have issues, but since this is a edge case prob best to go
        ## with the easiest solution (ie the latter)

        ## remap "to" intervals to reduced set
        to.r = reduce(to);
        if (length(to.r)<length(to)){
            to.map = values(gr.findoverlaps(to, to.r))[, c('query.id', 'subject.id')];
            values(to[to.map$query.id])$query.id = to.map$subject.id;            
        }

        # expand ranges as per mult
        rep.ix = as.integer(Rle(1:length(from), abs(values(from)$mult)))
        from = from[rep.ix]
        to = to[rep.ix]
        values(from)$mult = sign(values(from)$mult)
        
        ## determine seqlengths in the new genomeB
        tmp = aggregate(width(from), by = list(as.character(seqnames(to))), FUN = sum)
        new.lens = structure(seqlengths(genomeA), names = seqnames(genomeA))
        new.lens[tmp[,1]] = new.lens[tmp[,1]] + tmp[,2]                
        genomeB = Seqinfo(seqnames = names(new.lens), seqlengths = new.lens)
        
        ## determine loci of copied intervals in genomeB

        # to.tile is a tiling of genomeA including intervals that will and will not be mapped to genomeB
        # to.ix = NA means that interval WILL be mapped to genomeB
        to.tmp = gr.stripstrand(to); values(to.tmp)$to.ix = 1:length(to);
        to.gaps = gaps(to.tmp); to.gaps = to.gaps[strand(to.gaps)=='*']; values(to.gaps)$to.ix = NA;
        to.tile = sort(grbind(to.tmp, to.gaps));
        
        ## compute interval pairs by analyzing tile
        new.int = do.call('rbind', lapply(seqlevels(genomeA), function(x){
            this.tile = to.tile[seqnames(to.tile)==x]

            out = data.frame(chr.A = as.character(seqnames(this.tile)), start.A = start(this.tile),
              end.A = end(this.tile), strand.A = '+',
              chr.B = as.character(seqnames(this.tile)), start.B = start(this.tile),
              end.B = end(this.tile), strand.B = '+', flag = F, stringsAsFactors = F);

            insert.ix = !is.na(values(this.tile)$to.ix)
            if (any(insert.ix)){
                to.insert = to[values(this.tile)$to.ix[insert.ix]]
                from.insert = from[values(this.tile)$to.ix[insert.ix]]

                out$chr.A[insert.ix] = as.character(seqnames(from.insert))
                out$start.A[insert.ix] = start(from.insert)
                out$end.A[insert.ix] = end(from.insert)
                out$strand.A[insert.ix] = as.character(strand(from.insert));
                out$flag[insert.ix] = TRUE;
                
                width.A = out$end.A-out$start.A+1;
                out$start.B = cumsum(c(1,width.A[1:(length(width.A)-1)]))
                out$end.B = out$start.B + width.A - 1
                out$strand.B[insert.ix] = c('-', '+')[1+as.numeric(sign(values(from.insert)$mult)*c(-1, 1)[1+as.numeric(strand(to.insert)=='+')]>0)]
            }
            return(out)
        }));

        gc = gChain(GRanges(new.int$chr.A, IRanges(new.int$start.A, new.int$end.A), strand = new.int$strand.A, seqlengths = seqlengths(genomeA)),
            GRanges(new.int$chr.B, IRanges(new.int$start.B, new.int$end.B), strand = new.int$strand.B, seqlengths = seqlengths(genomeB)),
            val = data.frame(flag = new.int$flag))                  
        
    }                    

}


##########################################
# delete 
#
# creates a gChain that deletes the target interval(s) "target" (GRanges object
#
# the seqinfo for the deletion is taken from the "target" GRanges
##########################################
delete = function(target) {

    genomeA = seqinfo(target)
    old.lens = structure(seqlengths(genomeA), names = seqnames(genomeA))
      
    target = reduce(target);

    ## determine seqlengths in the new (reduced) genomeB
    tmp = aggregate(width(target), by = list(as.character(seqnames(target))), FUN = sum)
    new.lens = structure(seqlengths(genomeA), names = seqnames(genomeA))
    new.lens[tmp[,1]] = new.lens[tmp[,1]] - tmp[,2]                
    genomeB = Seqinfo(seqnames = names(new.lens), seqlengths = new.lens)
  
    ## to.tile is a tiling of genomeA including intervals that will and will not be mapped to genomeB
    ## to.ix = NA means that interval WILL be mapped to genomeB
    tgt.tmp = gr.stripstrand(target); 
    intA = sort(gaps(tgt.tmp), decreasing = F); intA = intA[strand(intA)=='*'];
  
    tmp = intA;
    tmp$id = 1:length(intA);
    tmp = sort(grbind(tmp, tgt.tmp))
    del.ix = which(is.na(tmp$id))
    flag.id = tmp$id[unique(c(pmax(1, del.ix-1), pmin(length(tmp), del.ix+1)))]
    flag = rep(F, length(intA))
    flag[flag.id] = TRUE

    strand(intA) = '+'
    intB = gr.fix(intA, genomeB);

    ix = unlist(split(1:length(intA), as.character(seqnames(intA))))
    starts = unlist(lapply(split(width(intA), as.character(seqnames(intA))), function(x){
        if (length(x)==1){
            1
        } else{
            cumsum(c(1,x[1:length(x)-1]))
        }
    }))
  
    starts[ix] = starts;
    ends = starts+width(intA)-1;

    intB = GRanges(seqnames(intA), IRanges(starts, ends), strand = '+', seqlengths = seqlengths(genomeB))

    return(gChain(intA, intB, val = data.frame(flag = flag)))

}
  



##########################################
# invert
#
# creates a gChain that inverts the target interval(s)
#
# if target intervals overlap --> then will flip each of the
# disjoint overlaps
#
# to do nested inversions .. use gChain composition of simple
# (nonoverlapping) inversions
#
##########################################
invert = function(target){

    genomeA = seqinfo(target)
    old.lens = structure(seqlengths(genomeA), names = seqnames(genomeA))
      
    tmp.target = disjoin(target);
    if (length(tmp.target)<length(target)){
        warning('Inverted intervals intersect .. disjoining and inverting individually')
    }
    target = tmp.target;

    ## determine seqlengths in the new (reduced) genomeB
    tmp = aggregate(width(target), by = list(as.character(seqnames(target))), FUN = sum)
    new.lens = structure(seqlengths(genomeA), names = seqnames(genomeA))
    new.lens[tmp[,1]] = new.lens[tmp[,1]] - tmp[,2]                
    genomeB = Seqinfo(seqnames = names(new.lens), seqlengths = new.lens)

    ## to.tile is a tiling of genomeA including intervals that will and will not be mapped to genomeB
    ## to.ix = NA means that interval WILL be mapped to genomeB
    tgt.tmp = gr.stripstrand(target);  values(tgt.tmp)$invert = T;
    intA = gaps(tgt.tmp); intA = intA[strand(intA)=='*']; values(intA)$invert = F;
    intA = sort(c(tgt.tmp, intA), decreasing = F);
    flag = width(tgt.tmp) < old.lens[as.character(seqnames(tgt.tmp))]
    strand(intA) = '+'
    intB = gr.fix(intA, genomeB)
    strand(intB)[values(intB)$invert] = '-'
  
    return(gChain(intA, intB, val = data.frame(flag = values(intB)$invert)))

}



###############################################
# permute 
#
# creates gChain to implement a genomic permutation, which is specified
# as a grl representing permutation cycles.  None of the intervals in the cycles
# can be overlapping each other.  
#
# ie each grl item i is a permutation cycle.  Every gr[j] in grl[[i]] item will be swapped with
# gr j+1  (or with gr[1] if j+1 > length(gr))
#
# Intervals in the genome that are not specified as part of a cycle will be preserved in their order
# (ie the same effect as if they were specified as part of a length 1 cycle)
###############################################
permute = function(cycles){

    if (!inherits(cycles, 'GRangesList')){
        cycles = GRangesList(cycles);
    }
        
    c.gr = gr.fix(unlist(cycles));
    values(c.gr)$ix = unlist(sapply(1:length(cycles), function(x) rep(x, length(cycles[[x]]))))

    if (any(strand(c.gr)=='*')){
        warning('Converting * strands to +')
        ix = which(as.logical(strand(c.gr)=='*'))
        strand(c.gr[ix]) = '+'
    }
    
    if (length(reduce(c.gr)) < length(c.gr)){
        stop('Error: Intervals describing permutation cycles have to be disjoint')
    }
    
    c.gr.gaps = gaps(c.gr, start = 1);
    c.gr.gaps = c.gr.gaps[which(as.logical(strand(c.gr.gaps)=='+'))]

    tile = sort(grbind(c.gr, c.gr.gaps))
    
    tile.ix = which(!is.na(tile$ix));
    cycles.ix = split(tile.ix, c.gr$ix[tile.ix])

    ## these are in the indices that we will replace with
    target.ix = levapply(tile.ix, c.gr$ix[tile.ix], function(x) x[((1:length(x)) %% length(x))+1])
    
    intA = tile;
    intA[tile.ix] = intA[target.ix];
    intB = gr.refactor(intA, seqnames(tile))
    strand(intB) = strand(tile);

    return(gChain(intA, intB, val = data.frame(flag = !is.na(intA$ix))))
}



################################################
# rearrange
#
# returns gChain applying a rearrangement to a genome specified as an "event", which is a GRanges of signed intervals.
# The event is specified as follows:
# The range to the right of a "+" strand event[i] will be connected to the left of a "+" strand event[i+1]
# (flip left / right for "-" strand in either of the above cases)
# If closed = T, then a connection will be made between event[length[event]] and event[1]
#
################################################
rearrange = function(event, ## this is a GRanges representing breakpoints comprising an event
    closed = TRUE, ## a "closed" event will connect the last breakpoint to the first
    retain = 0, ## whether or not to retain a breakpoint (i.e. an amplification bridge) as part of the left or right breakpoint
    filter.tel = TRUE, # these intervals (e.g. telomeres, centromeres) represent essential intervals required to keep the chromosome following the transformation
    filter.seg = NULL # if not null will also only filter contigs that map to a filter.seg segment (eg centromeres)  
)
{
    
    if (any(strand(event) == '*')){
        stop('bp1 and bp2 must be signed intervals (i.e. either + or -)')
    }
    
    if (sum(width(reduce(event))) != sum(width(event))){
        stop('Event cannot have duplicates with respect to location and strand')
    }
    
    values(event)$retain = retain    
    bp1 = event[-length(event)]
    bp2 = gr.flipstrand(event[-1])

    if (closed){
        bp1 = c(bp1, event[length(event)])
        bp2 = c(bp2, gr.flipstrand(event[1]))
    }
        
    sgn1 = c('-'=-1, '+'=1)[as.character(strand(bp1))]
    sgn2 = c('-'=-1, '+'=1)[as.character(strand(bp2))]

    event$is.bp = TRUE
    
    # first we tile the genome around the combined breakpoints
    g = gaps(gr.stripstrand(sort(event))); g = g[strand(g)=='*']; strand(g) = '+';
    g$is.bp = FALSE
    tile = grbind(event, g);
    tile = tile[order(gr.stripstrand(tile))]

    # label telomeric segments
    tile = gr.fix(tile);
    values(tile)$begin.tel = start(tile)==1
    values(tile)$end.tel = end(tile) == seqlengths(tile)[as.character(seqnames(tile))]
    values(tile)$is.tel = values(tile)$begin.tel | values(tile)$end.tel 
    values(tile)$tile.id = 1:length(tile);

    # label filter segments
    if (!is.null(filter.seg)){
        tile$filter.seg = seg.on.seg(tile, filter.seg)
    }      

    tmp = gr.findoverlaps(bp1, tile)
    bp1.ix = rep(NA, length(bp1))
    bp1.ix[tmp$query.id] = tmp$subject.id

    tmp = gr.findoverlaps(bp2, tile)
    bp2.ix = rep(NA, length(bp2))
    bp2.ix[tmp$query.id] = tmp$subject.id
    
    # now trim breakpoints by width(units) and add this
    # width to the interval pointed to by the breakpoint
    #
    # this will facilitate deletion and amplification "bridges"
    bpix = which(tile$is.bp)
    not.last = bpix<length(tile)
    not.first = bpix>1
    pos.bp = as.logical(strand(tile[bpix])== "+");
    wids = width(tile)
    wids[!values(tile)$retain] = 1;

    if (any(pos.bp)){
        if (any(end(tile)[bpix[pos.bp]]-wids[bpix[pos.bp]] > .Machine$integer.max)){
            warning('Trimmed some endpoints to below integer maximum')
        }
    }

    if (any(pos.bp & not.last)){
        if (any((start(tile)[bpix[pos.bp & not.last]+wids[bpix[pos.bp & not.last]]]-1) > .Machine$integer.max)){
            warning('Trimmed some endpoints to below integer maximum')
        }
    }
    
    if (any(!pos.bp)){
        if (any(start(tile)[bpix[!pos.bp]]+1 > .Machine$integer.max)){
           warning('Trimmed some endpoints to below integer maximum')
        }
    }

    if (any(!pos.bp & not.first)){
        if (any((end(tile)[bpix[!pos.bp & not.first]-wids[bpix[!pos.bp & not.first]]]+1) > .Machine$integer.max)){
           warning('Trimmed some endpoints to below integer maximum')
        }
    }
               
    end(tile)[bpix[pos.bp]] = pmax(start(tile)[bpix[pos.bp]], pmin(.Machine$integer.max, end(tile)[bpix[pos.bp]]-wids[bpix[pos.bp]]))
    start(tile)[bpix[pos.bp & not.last]+1] = pmin(.Machine$integer.max, start(tile)[bpix[pos.bp & not.last]+wids[bpix[pos.bp & not.last]]]-1)
    start(tile)[bpix[!pos.bp]] = pmin(.Machine$integer.max, start(tile)[bpix[!pos.bp]]+1);
    end(tile)[bpix[!pos.bp & not.first]-1] = pmin(.Machine$integer.max, end(tile)[bpix[!pos.bp & not.first]-wids[bpix[!pos.bp & not.first]]]+1);
    

    ab.pairs = cbind(pmax(1, pmin(length(tile), bp1.ix + sgn1)), pmax(0, pmin(length(tile), bp2.ix +sgn2)))
    
    pp = (sgn1*sgn2)>0 & sgn1>0;
    mm = (sgn1*sgn2)>0 & sgn1<0;
    mp = sgn1>0 & sgn2<0
    ab.pairs[pp,1] = -ab.pairs[pp,1] ## ++ breakpoints --> (-b)d adjacency
    ab.pairs[mm,2] = -ab.pairs[mm,2] ## -- breakpoints --> a(-c) adjacency
    ab.pairs[mp, ] = -ab.pairs[mp, ] ## +- breakpoints --> (-b)(-c) adjacency
    
    ## clean up adj pairs
    ## remove any that have crossed a chromosome boundary from their breakpoint
    ## this will occur in cases of badly formed breakpoint input (eg breakpoints that point outward
    ## from their telomeres)        
    keep = as.logical((seqnames(tile)[abs(ab.pairs[,1])]==seqnames(bp1)) & (seqnames(tile)[abs(ab.pairs[,2])]==seqnames(bp2))) & !tile$is.bp[abs(ab.pairs[,1])] & !tile$is.bp[abs(ab.pairs[,2])]
    ab.pairs = ab.pairs[keep, , drop = FALSE];
    ab.pairs = rbind(ab.pairs, cbind(-ab.pairs[,2], -ab.pairs[,1]));
    
    ## flag any segments flanking a breakpoint
    left.good = as.logical(seqnames(tile)[pmin(length(tile), bpix-1)] == seqnames(tile)[bpix])
    right.good = as.logical(seqnames(tile)[pmin(bpix+1, length(tile))] == seqnames(tile)[bpix])
    values(tile)$flag = (1:length(tile)) %in% c(bpix[left.good]-1, bpix[right.good]+1)
    
    ## build "aberrant" adjacency matrix representing directed graph of edges connecting
    ## <signed> nodes.
    adj.ab = matrix(FALSE, nrow = 2*length(tile), ncol = 2*length(tile), dimnames = rep(list(as.character(c(1:length(tile), -(1:length(tile))))), 2)) 
    tmp = table(ab.pairs[,1], ab.pairs[,2])

    if (nrow(tmp)>0){
        ix = which(tmp!=0, arr.ind = TRUE);
        ix.n = cbind(rownames(tmp)[ix[,1]], colnames(tmp)[ix[,2]])
        adj.ab[ix.n] = tmp[ix.n];
    }
    
    ## build reference adjacency matrix (representing consecutive segments on the reference genome)
    seg.ix = which(!tile$is.bp)
    ref.pairs = cbind(seg.ix[1:(length(seg.ix)-1)], seg.ix[2:(length(seg.ix))])
    ref.pairs = ref.pairs[ref.pairs[,1]>0 & ref.pairs[,2]!=length(tile), , drop = FALSE]
    ref.pairs = ref.pairs[which(as.logical(seqnames(tile[ref.pairs[,1]]) == seqnames(tile[ref.pairs[,2]]))), , drop = FALSE]
    ref.pairs = rbind(ref.pairs, cbind(-ref.pairs[,2], -ref.pairs[,1])) # reverse ref pairs
      
    adj.ref = matrix(FALSE, nrow = 2*length(tile), ncol = 2*length(tile), dimnames = rep(list(as.character(c(1:length(tile), -(1:length(tile))))), 2))
    tmp = table(ref.pairs[,1], ref.pairs[,2])
    ix = which(tmp!=0, arr.ind = TRUE);
    ix.n = cbind(rownames(tmp)[ix[,1]], colnames(tmp)[ix[,2]])
    adj.ref[ix.n] = tmp[ix.n];

    # traverse all non-cyclic paths through adj.ab (ie no nodes repeated) seeded by source nodes in adj.ref
    # (ie telomere ends)
    #
    # adj.ref are used ONLY when there is no non-aberrant connection in a segment
    # (note - this will only happen at "internal" segments) .. since + q telomeres and - p telomeres will
    # have no outgoing edges)    
    # not most efficient impl but should do the trck (FIX: consider using igraph when back up and running)
    seeds = which(colSums(adj.ref)==0 & !values(tile)$is.bp[abs(as.numeric((rownames(adj.ref))))])
    I = diag(rep(TRUE, nrow(adj.ab)))
    paths = copies = list();
    for (i in seeds){
        these.paths = list(i)
        these.copies = list(1+adj.ab[i,i])
        tail = i; # this is vector containing tail item of these.paths
        visited = array(FALSE, dim = c(1, ncol(adj.ab)))   ## visited is paths x nodes matrix keeps track of visited nodes in each path
        done = FALSE ## paths are done if their "tail" node has no non-self children that have not already been visited
        ## cycle through !done paths
        while (any(!done)){
            j = which(!done)[1];
            i = tail[j]
            visited[j, i] = TRUE;

            children = which(adj.ab[i, , drop = FALSE]>0 & I[i,, drop = FALSE]==0 & !visited[j, ])
            

            if (length(children)>0){
                these.paths = c(these.paths[-j], lapply(children, function(x) c(these.paths[[j]], x)))
                these.copies = c(these.copies[-j], lapply(children, function(x) c(these.copies[[j]], 1+adj.ab[x,x])));
                visited = rbind(visited[-j, , drop = FALSE], visited[rep(j, length(children)), , drop = FALSE])
                tail = c(tail[-j], children);
                done = c(done[-j], done[rep(j, length(children))]);                
            } else{
                done[j] = TRUE
            }
        }        

        paths = c(paths, these.paths);
        copies = c(copies, these.copies);

    }
    
    # contig data structure will have new chromosomes and orders of segments on chromosomes    
    contigs = data.frame(node.id = rownames(adj.ab)[unlist(paths)], contig.id = unlist(lapply(1:length(paths), function(x) rep(x, length(paths[[x]])))), ord = unlist(lapply(paths, function(x) 1:length(x))), copies = unlist(copies), tile.id = abs(as.numeric(rownames(adj.ab)[unlist(paths)])), sign = sign(as.numeric(rownames(adj.ab)[unlist(paths)])), stringsAsFactors = F)
            
    # now dedup by identifying contigs that are reverse complements of each other
    ctags = sapply(split(contigs[, 'node.id'], contigs[, 'contig.id']), function(x) paste(x, collapse = ', '))
    ctags.rev = sapply(split(contigs[, 'node.id'], contigs[, 'contig.id']), function(x) paste(rev(-as.numeric(x)), collapse = ', '))
    m = match(ctags, ctags.rev)
    keep = sapply(1:length(m), function(x) if (is.na(m[x])) T else m[x]>=x)
    contigs = contigs[contigs$contig.id %in% as.numeric(names(ctags)[which(keep)]), ]

    ## keep only contigs that begin and end in a tel
    if (filter.tel) {
        has.begin = vaggregate(as.numeric(node.id) ~ contig.id, data = contigs, function(x) if (x[1]>0) tile$begin.tel[abs(x[1])] else tile$end.tel[abs(x[1])])
        has.end = vaggregate(as.numeric(node.id) ~ contig.id, data = contigs, function(x) if (x[length(x)]>0) tile$end.tel[abs(x[length(x)])] else tile$begin.tel[abs(x[length(x)])])
        contigs = contigs[contigs$contig.id %in% intersect(names(which(has.begin)), names(which(has.end))), ]
    }

    if (!is.null(filter.seg)){
        good.contigs = unique(contigs$contig.id[contigs$tile.id %in% tile$tile.id[tile$filter.seg]]);
        contigs = contigs[contigs$contig.id %in% good.contigs, ]
    }

    cnames.og = structure(as.character(seqnames(tile)[contigs$tile.id[contigs$ord == 1]]),
        names = contigs$contig.id[contigs$ord == 1])
    cnames = structure(dedup(as.character(seqnames(tile)[contigs$tile.id[contigs$ord == 1]]), ' der '),
        names = contigs$contig.id[contigs$ord == 1])
    contigs$cname = cnames[as.character(contigs$contig.id)]
    contigs$cname.og = cnames.og[as.character(contigs$contig.id)]
    contigs = contigs[order(match(contigs$cname.og, seqlevels(bp1))), ];
    
    # now make maps
    ix.expand = as.numeric(Rle(1:nrow(contigs), contigs$copies))
    intA = tile[contigs$tile.id[ix.expand]]
    intB = gr.refactor(intA, contigs$cname[ix.expand])
    strand(intB)[contigs$sign[ix.expand]<0] = '-'    

    return(gChain(intA, intB, val = as.data.frame(values(intA)[, 'flag', drop = FALSE])))
}





###########################
# bfb
#
# makes gChain simulating a random BFB on chromosome "chr" with bridges occuring
# within "w" fraction (if 0 < w < 1) or bases (if w > 1) of end of "chr" (if end = T, start otherwise)
# and new breakages occuring within "wd" fraction (if 0 < wd < 1) or "wd" bases (if wd > 1) of chromosome end
#
#
###########################
bfb = function(chr, numcycles = 5, w = 0.2, wd = 0.4, end = TRUE, si = gUtils::si, verbose = TRUE, max.chrom = 800e6){
    if (w<0 | wd < 0){
        stop('w and wd must be positive')
    }
          
    gc = gChain(seqinfo2gr(si))
    si1 = seqinfo2gr(seqinfo(gc)[[1]])
    
    for (i in 1:numcycles){
        if (verbose){
            message('Cycle ', i, ' chrom width ', seqlengths(seqinfo(gc)[[2]])[chr]/1e6, ' MB\n')
        }
        si1 = seqinfo2gr(seqinfo(gc)[[2]])
        si.this = si1[chr]
        si1.other = si1[setdiff(seqlevels(si1), chr)]        
        
        if (w<1){
            W = width(si.this)*w
        } else{
            W = pmin(w, width(si.this))
        }
        
        if (end){
            bp = pmax(1, width(si.this)-runif(1)*W)
            inv = GRanges(chr, IRanges(1, bp), strand = '+', seqlengths = seqlengths(si1))
            tmp.gr = c(inv, gr.flipstrand(inv), si1.other);
            gc1 = gChain(c(inv, inv, si1.other), gr.refactor(tmp.gr, seqnames(tmp.gr)), val = data.frame(flag = 1:length(tmp.gr) %in% 1:2))
        } else{
            bp = runif(1)*W
            inv = GRanges(chr, IRanges(bp, width(si.this)), strand = '+', seqlengths = seqlengths(si1))
            tmp.gr = c(gr.flipstrand(inv), inv, si1.other);
            gc1 = gChain(c(inv, inv, si1.other), gr.refactor(tmp.gr, seqnames(tmp.gr)), val = data.frame(flag = 1:length(tmp.gr) %in% 1:2))
        }  
                    
        si2 = seqinfo2gr(seqinfo(gc1)[[2]])
        
        if (wd<1){
            WD = width(si2[chr])*wd
        } else{
            WD = wd
        }
          
        if (end){
            to.del = gr.end(si2[chr], pmax(seqlengths(si2)[chr] - max.chrom, round(runif(1)*WD)))
        } else{
            to.del = gr.start(si2[chr], pmax(seqlengths(si2)[chr] - max.chrom, round(runif(1)*WD)))
        }
        
        gc2 = delete(to.del)
        
        gc = gc2*gc1*gc
    }

    return(gc)
}

########################
# utility functions for gChain
#
#
########################
          
## vec2ir
##
## given (integer) vector v converts sequential runs of 1 in diff(x) into
## ranges describing start and end positions of those runs.
##
## eg a vector x = c(1,2,3,10,11,12,5,6,7) would be summarized as IRanges(start = c(1, 10, 5), end = c(3, 12, 7))
##
## a vector x = seq(1, 100, 3) (ie with no sequential runs) would be described as IRanges(start = x, end = x)
##
vec2ir = function(x){
    bkpoints.x = which(diff(floor(x))!=1)
    return(IRanges(start = x[c(1, bkpoints.x+1)], end = x[c(bkpoints.x, length(x))]))
}

## ir2vec
##
## returns vector "expanding" the positions specified by the given IRanges input x
## optional scalar or vector arg "rev" will optionally reverse sequences for the ranges
## for which rev is TRUE
##
## if "each" is scalar or vector of length(x) will replicate each item of the output vector resulting from x[k]
## "each" or "each[k]" times, respctively
ir2vec = function(x, rev = FALSE, each = NULL){
    if (length(rev) == 1){
        rev = rep(rev, length(x))
    }
    if (!is.null(each)){
        if (length(each) == 1){
            each = rep(each, length(x))
        }
        return(unlist(lapply(1:length(x),
            function(i) if (rev[i]) rep(rev(c(start(x[i]):end(x[i]))), each = each[i])
                else rep(c(start(x[i]):end(x[i])), each = each[i]))))
    } else{
        return(unlist(lapply(1:length(x), function(i) if (rev[i]) rev(c(start(x[i]):end(x[i]))) else c(start(x[i]):end(x[i])))))
    }   
}



### concatentate gChains
gCat <- function(x, ...) {

    if (missing('x')){
        args <- list(...)
    } else{
        args <- c(x, list(...))
    }

    ## remove empty
    args <- args[sapply(args, function(y) length(y@.galx)) > 0]

    print('gCat: working on data tables')
    dts <- lapply(args, function(x) {
        sn <- as.character(seqnames(x@.galx))
        if (identical('NA', sn)){
            return(data.table())
        }
        st <- start(x@.galx)
        ed <- end(x@.galx)
        sr <- as.character(strand(x@.galx))
        dt <- data.table(seqnames=sn, start=st, end=ed, strand=sr)
        return(dt)
    })

    dtx <- rbindlist(dts)
  
    dts <- lapply(args, function(x) {
        sn <- as.character(seqnames(x@.galy))
        if (identical('NA', sn)){
            return(data.table())
        }
        st <- start(x@.galy)
        ed <- end(x@.galy)
        sr <- as.character(strand(x@.galy))
        dt <- data.table(seqnames=sn, start=st, end=ed, strand=sr)
        return(dt)
    })

    dty <- rbindlist(dts)

    slx <- do.call('c', lapply(args, function(y) seqlengths(y@.galx)))
    slx <- slx[unique(names(slx))]
    sly <- do.call('c', lapply(args, function(y) seqlengths(y@.galy)))
    sly <- sly[unique(names(sly))]

    if (nrow(dtx) == 0){
        return(gChain())
    }

    ## I hacked the GRanges constructor to make it way faster for large dt to Gr conversions
    print('...gCat: making GRanges')
    dtx[, len := max(end), by='seqnames']
    gr.dtx <- GRanges()
    gr.dtx@ranges <- IRanges(dtx$start, dtx$end)
    gr.dtx@seqnames <- Rle(factor(dtx$seqnames, levels=unique(dtx$seqnames)))
    gr.dtx@strand <- Rle(factor(dtx$strand, levels=c('+', '-', '*')))
    df <- DataFrame(rep(1, nrow(dtx)))
    gr.dtx@elementMetadata <- df[,c()]
    gr.dtx@seqinfo <- Seqinfo(as.character(unique(dtx$seqnames)), seqlengths=dtx$len[!duplicated(dtx$seqnames)])
  
    gr.dty <- GRanges()
    dty[, len := max(end), by='seqnames']  
    gr.dty@ranges <- IRanges(dty$start, dty$end)
    gr.dty@seqnames <- Rle(factor(dty$seqnames, levels=unique(dty$seqnames)))
    gr.dty@strand <- Rle(factor(dty$strand, levels=c('+', '-', '*')))
    df <- DataFrame(rep(1, nrow(dty)))
    gr.dty@elementMetadata <- df[,c()]
    gr.dty@seqinfo <- Seqinfo(as.character(unique(dty$seqnames)), seqlengths=dty$end[!duplicated(dty$seqnames)])  

    pad.left  <- unlist(lapply(args, function(x) x@.pad.left))
    pad.right <- unlist(lapply(args, function(x) x@.pad.right))

    print('...gCat: making gChain')
    return(gChain(gr.dtx, gr.dty, pad.left=pad.left, pad.right=pad.right))
  
}



#' @name gUnique
#' @title gUnique
#' @description
#'
#' remove duplicate x to y mappings
#'
#' Will select only the longest, and most "leftward" in genomic coordaintes
#' (in that order) entry of each unique x to y pair, where
#' uniqueness is determined by seqnames and strand, but not by position.
#' This could occur, for instance, if a read has multiple mappings to a
#' genomic contig or multiple mappings within one chromosome. Occasionally,
#' one might want to take only the longest match, and discard secondary alignments.
#'
#' @param gc gChain object to deduplicate
#' @return Deduplicated gChain
#' @note This is a bit of a hack, as there are better ways to remove dupes before this step
#' @export
gUnique = function(gc){

    lix <- links(gc)$x
    liy <- links(gc)$y
    snx <- as.character(seqnames(lix))
    sny <- as.character(seqnames(liy))
    stx <- as.character(strand(lix))
    sty <- as.character(strand(liy))
    wid <- width(lix)
    start <- start(lix)

    ## Need to do -wid so that longest is first
    gcx <- data.table(id=seq_along(lix), snx=snx, sny=sny, stx=stx, sty=sty, wid=-wid, start=start)
    setkey(gcx, wid, start, snx, sny, stx, sty) # sort by width, then start
    setkey(gcx, snx, sny, stx, sty) # remove wid and start as keys
    gcx <- unique(gcx)

    suppressWarnings(out <- new('gChain', x=lix[gcx$id], y=liy[gcx$id]))
  
    return(out)
}



#' @name squeeze
#' @title squeeze
#' @description
#'
#' "squeezes" pile of IRanges so that width(squeeze(ir)) = width(ir), end(ir[length(ir)]) = sum(width(ir)) + 1
#' start(squeeze(ir)[k]) = end(squeeze(ir)[k-1])+1 for all k>1, and start(squeeze(ir))[1] = 1
#'
#' optional extra input "gap" yields following output: width(squeeze(ir, gap)) = width(ir), end(ir[length(ir)]) = sum(width(ir)) + gap*(length(ir)-1) + 1
#' start(squeeze(ir, gap)[k]) = end(squeeze(ir, gap)[k-1])+1+gap for all k>1, and start(squeeze(ir, gap))[1] = 1
#'
#' @param x IRanges object
#' @param gap integer bp gap
#' @return squeezed IRanges object
squeeze = function(x, gap = 0){

    if (!inherits(x, 'IRanges')){
        stop('squeeze() only defined for IRanges')            
    }
    if (length(x) == 0){
        return(x)
    } else if (length(x) == 1){
        return(IRanges(start = 1, width = width(x)))
    } else{
        starts = cumsum(c(1, width(x[1:(length(x)-1)])+gap))
        ends = starts+width(x)-1
        return(IRanges(starts, ends))
    }
}




#################################
#' @name gSubset
#' @title gSubset
#' @description
#' Subset a gChain
#'
#' Make smaller gChain that contains only the seqnames
#' given by xnames and/or ynames.
#' @param gc gChain to subset
#' @param xnames character vector that has seqnames of links(gc)$x to keep
#' @param ynames character vector that has seqnames of links(gc)$y to keep
#' @param x.or.y Links with hits in xnames OR ynames included. Default FALSE (i.e. AND)
#' @return subsetted gChain
#' @export
gSubset = function(gc, xnames=NULL, ynames=NULL, x.or.y = FALSE){

    if (class(gc) != 'gChain'){
        stop('gSubset: Need to input a gChain')
    }

    # default is to include everything
    val.x <- val.y <- rep(TRUE, length(links(gc)$x))
    if (!is.null(xnames)){
        val.x <- as.character(seqnames(links(gc)$x)) %in% xnames
    }
    if (!is.null(ynames)){
        val.y <- as.character(seqnames(links(gc)$y)) %in% ynames
    }
    if (x.or.y){
        val <- val.x | val.y
    } else{
        val <- val.x & val.y
    } 

    sn <- as.character(seqnames(links(gc)$y))[val]
    start <- start(links(gc)$y)[val]
    end <- end(links(gc)$y)[val]
    strand <- as.character(strand(links(gc)$y))[val]
    gc@.galy <- GRanges(sn, IRanges(start, end), strand=strand)

    sn <- as.character(seqnames(links(gc)$x))[val]
    start <- start(links(gc)$x)[val]
    end <- end(links(gc)$x)[val]
    strand <- as.character(strand(links(gc)$x))[val]
    gc@.galx <- GRanges(sn, IRanges(start, end), strand=strand)
  

    gc@.pad.left  <- gc@.pad.left[val] 
    gc@.pad.right <- gc@.pad.right[val]
    gc@.scale <- gc@.scale[val]
    if (nrow(gc@values) != 0){
        gc@values <- gc@values[val, ]
    }
    gc@.n <- sum(width(gc@.galx))
    gc@.m <- sum(width(gc@.galy))  

    return(gc)
}




#' @name grl.split
#' @title grl.split
#' @description
#'
#' splits GRL's with respect to their seqnames and strand (default), returning
#' new grl whose items only contain ranges with a single seqname / strand
#'
#' can also split by arbitrary (specified) genomic ranges value fields
#'
#' @param grl \code{GRangesList} to split
#' @param seqname Default TRUE
#' @param strand Default TRUE
#' @param values columns of values field in grl
grl.split = function(grl, seqname = TRUE, strand = TRUE, values = c())
{
    ele = tryCatch(as.data.frame(grl)$element, error = function(e) e)
    if (inherits(ele, 'error') | is.null(ele)){
        if (is.null(names(grl))){
            nm = 1:length(names(grl))
        } else{
            nm = names(grl)
        }

        ele = unlist(lapply(1:length(grl), function(x) rep(nm[x], length(grl[[x]]))))
    }

    gr = unlist(grl)
    names(gr) = NULL;

    by = ele;
    if (seqname){
        by = paste(by, seqnames(gr))
    }

    if (strand){
        by = paste(by, strand(gr))
    }

    values = intersect(names(values(gr)), values);

    if (length(values)>0){
        for (val in values){
            by = paste(by, values(gr)[, val])
        }
    }

    out = split(gr, by);
    names(out) = ele[!duplicated(by)]

    values(out) = values(grl[ele[!duplicated(by)]])

    return(out)
}


vaggregate = function(...){
    out = aggregate(...);
    return(structure(out[,ncol(out)], names = do.call(paste, lapply(names(out)[1:(ncol(out)-1)], function(x) out[,x]))))
}


levapply = function(x, by, FUN = 'order'){
    if (!is.list(by)){
      by = list(by)
    }

    f = factor(do.call('paste', c(list(sep = '|'), by)))
    ixl = split(1:length(x), f);
    ixv = lapply(ixl, function(y) x[y])
    res = structure(unlist(lapply(ixv, FUN)), names = unlist(ixl))
    out = rep(NA, length(x))
    out[as.numeric(names(res))] = res;
    return(out)
}



seqinfo2gr = function(si, strip.empty = FALSE){
    ## treat si as seqlengths if vector
    if (is(si, 'vector')){
        si = Seqinfo(seqlengths = si, seqnames = names(si))
    } else if (!is(si, 'Seqinfo')){
        si = seqinfo(si)
    }

    sl = seqlengths(si)
    sn = seqnames(si);
    sl[is.na(sl)] = 0;

    if (strip.empty){
        sn = sn[sl!=0];
        sl = sl[sl!=0];
    }

    sigr = GRanges(sn, IRanges(rep(1, length(sl)), width = sl), seqlengths = seqlengths(si), strand = rep('+', length(sl)))
    names(sigr) = sn;

    return(sigr)
}




dedup = function(x, suffix = '.'){
    dup = duplicated(x);
    udup = setdiff(unique(x[dup]), NA)
    udup.ix = lapply(udup, function(y) which(x==y))
    udup.suffices = lapply(udup.ix, function(y) c('', paste(suffix, 2:length(y), sep = '')))
    out = x;
    out[unlist(udup.ix)] = paste(out[unlist(udup.ix)], unlist(udup.suffices), sep = '');
    return(out)  
}





gr.pad = function(gr, pad){
    start(gr) = pmax(1, start(gr)-pad)
    en = pmin(seqlengths(gr)[as.character(seqnames(gr))], end(gr)+pad)
    end(gr) = ifelse(is.na(en), end(gr)+pad, en)
    return(gr)
}



#' gr.refactor
#'
#' Takes a pile of ranges gr and new seqnames "sn" (either of length 1 or
#' of length(gr)) and returns a gr object with the new seqnames and same
#' widths and new start coordinates.  These coordinates are determined by placing
#' each gr on the corresponding chromosome at 1 + gap after previous gr (or at 1)
#' @param gr \code{GRanges} to refactor
#' @param sn character vector of new seqnames
#' @param gap Default 0
#' @param rev Default FALSE
gr.refactor = function(gr, sn, gap = 0, rev = FALSE){
    
    if (is.factor(sn)){
        slev = levels(sn)
    }
    else{
        slev = unique(sn);
    }

    sn = cbind(as.character(start(gr)), as.character(sn))[,2]
    w = width(gr)
    gap = pmax(cbind(gap, w)[,1], 0);

    starts = levapply(1:length(w), sn, function(x) cumsum(gap[x] + c(1, w[x[1:(length(x)-1)]])[1:length(x)])-gap[x])
    ir = IRanges(starts, width = width(gr))

    sl = aggregate(end(ir)+gap, by = list(sn), FUN = max); sl = structure(sl[,2], names = sl[,1])

    oth.names = setdiff(slev, names(sl))
    if (length(oth.names)>0){
        sl[oth.names] = NA
    }
    sl = sl[slev]

    out = GRanges(sn, ir, strand = strand(gr), seqlengths = sl)
    values(out) = values(gr);

    return(out)
}


#' gr.tostring
#'
#' dumps out a quick text representation of a gr object (ie a character vector)
#'
#' @param gr \code{GRanges}
#' @param places Number of decimal places. Default 2
#' @param interval Default 1e6
#' @param unit Default "MB"
#' @param prefix Default "chr"
#' @return text representation of input
gr.tostring = function(gr, places = 2, interval = 1e6, unit = 'MB', prefix = 'chr'){
    p1 = round(start(gr)/interval, places);
    p2 = round(end(gr)/interval, places);
    return(paste(prefix, as.character(seqnames(gr)), ':', p1, '-', p2, ' ', unit, sep = ''));
}




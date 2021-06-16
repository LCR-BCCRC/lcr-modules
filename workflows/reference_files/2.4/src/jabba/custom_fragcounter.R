#' @import GenomicRanges
#' @import gUtils
#' @import rtracklayer
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom data.table data.table fread rbindlist set setkey setkeyv setnames transpose as.data.table
#' @importFrom stats cor loess predict quantile
#' @importFrom skidb read_gencode
#' @importFrom Biostrings alphabetFrequency
#' @importFrom parallel mclapply
#' @importFrom Rsamtools BamFile
#' @importFrom grDevices col2rgb rgb png dev.off
#' @importFrom IRanges IRanges
#' @importFrom DNAcopy CNA segment smooth.CNA
#' @importFrom graphics plot lines
#' @importFrom utils globalVariables



globalVariables(c(".", ".N", ".SD", ":=", "V1", "V2", "V3", "bam", "bin", "black_list_pct", "blacklisted", "child", "chr", "count", "counts", "gc.wig", "index", "ix.start", "lev0", "lev1", "map.wig", "newcount", "numlevs", "pair", "parent", "reads", "reads.corrected", "rowid", "sortSeqlevels", "tmp"))




#' @title Multi-scale coverage correction
#' @description Given gc and mappability coverage correction at k "nested" scales finds the coverage
#' assignment at the finest scale that yields the best correction at every scale
#' @name multicoco
#' @param cov GRanges constant with GRanges of coverage samples with (by default) fields $reads, $map, $gc
#' @param numlevs integer numbers of scales at which to correct
#' @param base integer Scale multiplier
#' @param fields character vector fields of gc to use as covariates
#' @param iterative boolean whether to iterate
#' @param presegment boolean whether to presegment
#' @param min.segwidth integer when presegmenting, minimum segment width
#' @param mono boolean Wether to only do single iteration at finest scale
#' @param verbose boolean Wether to print log to console
#' @param FUN function with which to correct coverage (by default loess
#' correction modified from HMMcopy that takes in granges with fields
#' $reads and other fields specified in "fields"
#' @param ... additional args to FUN
#' @param mc.cores integer Number of cores to use
#' @param exome boolean If TRUE, perform correction using exons as bins instead of fixed size
#' @author Marcin Imielinski
#' @usage multicoco(cov, numlevs = 1, base = max(10, 1e5/max(width(cov))),
#' fields = c("gc", "map"), iterative = TRUE, presegment = TRUE,
#' min.segwidth = 5e6, mono = TRUE, verbose = TRUE, FUN = NULL, ...,
#' mc.cores = 1, exome = FALSE)
#' @export


multicoco = function(cov, numlevs = 1, base = max(10, 1e5 / max(width(cov))),
                     fields = c("gc", "map"), iterative = TRUE, presegment = TRUE,
                     min.segwidth = 5e6, mono = TRUE, verbose = TRUE,
                     FUN = NULL, ..., mc.cores = 1, exome = FALSE) {
  if (verbose) {
    cat('Converting to data.table\n')
  }
  WID = max(width(cov))
##  library(data.table) #' twalradt Wednesday, Jan 16, 2019 03:54:24 PM
  cov.dt = gr2dt(cov)        
  sl = structure(as.numeric(1:length(GenomeInfoDb::seqlevels(cov))), names = GenomeInfoDb::seqlevels(cov))       
  if (verbose) {
    cat('Grouping intervals\n')
  }
  ## compute level means
  ## lev 0 = raw data
  ## lev 1 = base-fold collapsed
  ## lev 2 = base^2-fold collapsed
  ## etc
  parentmap= list() ## data.tables that map lev labels at level k  to parent lev labels
  cov.dt[, lev0:=as.character(1:length(seqnames))]
  for (k in 1:numlevs) {
    if (verbose) {
      cat('Defining', base^k, 'fold collapsed ranges\n')
    }
        
    ## cov.dt[, eval(paste("lev", k, sep = '')) := as.character(sl[seqnames] + as.numeric(Rle(as.character(1:length(start)), rep(base^k, length(start)))[1:length(start)])/length(start)), by = seqnames]


    ##  -.01 just to make Trent Happy ie produce identical result to above
    ## This line creates column that you can use to collapse the data with
    cov.dt[, eval(paste("lev", k, sep = '')) :=  as.character(sl[seqnames] + ceiling((1:.N-0.01)/base^k)), by = seqnames]


    parentmap[[k]] = data.table(parent = cov.dt[, get(paste("lev", k, sep = ''))], child = cov.dt[, get(paste("lev", k-1, sep = ''))], key = 'child')[!duplicated(child), ]
  }
  if (presegment) { ## perform rough segmentation at highest level
    seg = NULL
    sl = GenomeInfoDb::seqlengths(cov)
    if (verbose) {
      cat('Presegmenting at ', as.integer(WID*base^(numlevs)), ' bp scale \n')
    }
##    require(DNAcopy) #' twalradt Wednesday, Jan 16, 2019 04:36:46 PM
    set.seed(42) #' twalradt Friday, Apr 20, 2018 01:07:28 PM
    if (exome == TRUE) {
      tiles = gr.tile(sl, width = 1e5)
      tiles = tiles %Q% (seqnames %in% c(seq(1,22),"X","Y"))
      GenomeInfoDb::seqlevels(tiles) = GenomeInfoDb::seqlevelsInUse(tiles) # Get rid of empty seqname factor levels
      tmp.cov = tiles %$% cov
      tmp.cov = tmp.cov %Q% (!is.na(reads))
      # tmp.cov = gr.val(tiles, cov, val = c('reads', 'gc', 'map'), na.rm = TRUE)
    } else {
      tmp.cov = seg2gr(cov.dt[,list(chr = seqnames[1], start = min(start), end = max(end), strand = strand[1], reads = mean(reads, na.rm = T)), by = get(paste("lev", numlevs, sep = ''))][end>start, ], seqlengths = sl)
    }
    ix = which(!is.na(values(tmp.cov)[, 'reads']))
    tmp = data.frame()
    tryCatch({
      cna = CNA(log(values(tmp.cov)[, 'reads'])[ix], as.character(seqnames(tmp.cov))[ix], start(tmp.cov)[ix], data.type = 'logratio')
      tmp = print(segment(smooth.CNA(cna), alpha = 1e-5, verbose = T))
      tmp = tmp[!is.na(tmp$loc.start) & !is.na(tmp$chrom) & !is.na(tmp$loc.end), , drop = F]
    }, error = function(e) warning('DNACopy error moving on without segmenting'))
    if (nrow(tmp)>0) {
      seg = sort(seg2gr(tmp, seqlengths = sl))
      seg = seg[width(seg)>min.segwidth] ## remove small segments
      seg.dt = gr2dt(seg)
      if (nrow(seg.dt)>0) {
        seg = seg2gr(seg.dt[, list(seqnames = seqnames,
                                   start = ifelse(c(FALSE, seqnames[-length(seqnames)]==seqnames[-1]), c(1, start[-1]), 1),
                                   end = ifelse(c(seqnames[-length(seqnames)]==seqnames[-1], FALSE), c(start[-1]-1, Inf), GenomeInfoDb::seqlengths(seg)[as.character(seqnames)]))], seqlengths = sl)
        seg = gr.val(seg, tmp.cov, 'reads') ## populate with mean coverage
        seg$reads = seg$reads/sum(as.numeric(seg$reads*width(seg))/sum(as.numeric(width(seg)))) ## normalize segs by weigthed mean (so these become a correction factor)
      }
      else {
        seg = NULL
      }
    }
    else {
      seg = NULL
    }
  }
  else {
    seg = NULL
  }
  if (verbose) {
    cat('Aggregating coverage within levels \n')
  }
  ## list of data frames showing scales of increasing collapse
  cov.dt[, ix := 1:nrow(cov.dt)]
  cmd1 = paste('list(ix.start = ix[1], ix.end = ix[length(ix)], reads = mean(reads, na.rm = T),', paste(sapply(fields, function(f) sprintf("%s = mean(%s, na.rm = T)", f, f)), collapse = ','), ')', sep = '')
  cmd2 = paste('list(lab = lev0, reads,', paste(fields, collapse = ','), ', seqnames, start, end)', sep = '')
  if (mono) {
    if (verbose) {
      cat('Mono scale correction \n')
    }
    grs = list(cov.dt[, eval(parse(text=cmd2))])
    numlevs = 1
  }
  else {
    grs = c( list(cov.dt[, eval(parse(text=cmd2))]), lapply(1:numlevs, function(x) {
      if (verbose) {
        cat('Aggregating coverage in level', x,  '\n')
      }
      out = cov.dt[, eval(parse(text=cmd1)), keyby = list(lab = get(paste('lev', x, sep = '')))]
      out[, ":="(seqnames = cov.dt$seqnames[ix.start], end = cov.dt$end[ix.start], start = cov.dt$start[ix.start])]
      out[, ":="(ix.start= NULL, ix.end = NULL)]
      return(out)
    }))
  }
  setkey(grs[[1]], 'lab')
  ## modified from HMMCopy to
  ## (1) take arbitrary set of covariates, specified by fields vector
  ## (2) employ as input an optional preliminary (coarse) segmentation with which to adjust signal immediately prior to loess
  ## NOTE: (this only impacts the loess fitting, does not impose any segmentation on the data)
  ##
  if (is.null(FUN)) {
    FUN = function(x, fields = fields, samplesize = 5e4, seg = NULL, verbose = T, doutlier = 0.001, routlier = 0.01) {
      if (!all(fields %in% names(x))) {
        stop(paste('Missing columns:', paste(fields[!(fields %in% names(x))], collapse = ',')))
      }
      x$valid <- TRUE
      for (f in fields) {
        x$valid[is.na(x[, f])] = FALSE
        x$valid[which(is.infinite(x[, f]))] = FALSE
      }
      if (verbose) {
        cat('Quantile filtering response and covariates\n')
      }
      range <- quantile(x$reads[x$valid], prob = c(routlier, 1 - routlier), na.rm = TRUE)
      if (verbose) {
        cat(sprintf("Response min quantile: %s max quantile: %s\n", round(range[1],2), round(range[2],2)))
      }
      domains = lapply(fields, function(f) quantile(x[x$valid, f], prob = c(doutlier, 1 - doutlier), na.rm = TRUE))
      names(domains) = fields
      x$ideal <- x$valid
      x$ideal[x$reads<=range[1] | x$reads>range[2]] = FALSE
      for (f in fields)
      {
        x$ideal[x[, f] < domains[[f]][1] | x[, f] > domains[[f]][2]] = FALSE
      }
      if (verbose) {
        cat(sprintf('Nominated %s of %s data points for loess fitting\n', sum(x$ideal), nrow(x)))
      }
      set <- which(x$ideal)
      if (length(set)<=10) {
        warning("Not enough samples for loess fitting - check to see if missing or truncated data?")
        return(x$reads)
      }
      for (f in fields) {
        if (verbose) {
          message(sprintf("Correcting for %s bias...", f))
        }
        set.seed(42)
        select <- sample(set, min(length(set), samplesize))
        x2 = x[, c('reads', f)]
        x2$covariate = x2[, f]
        x2s = x2[select, ]                    
        if (!is.null(seg)) {  ## here we apply a prelmiinary segmentation to correct for large scale copy variation allow more power to reveal the covariate signal
          if (verbose) {
            message('Applying preliminary segmentation prior to loess fitting')
          }
          x.grs = gr.val(seg2gr(x[select, ], seqlengths = NULL), seg, 'reads')
          x2s$reads = x2s$reads/x.grs$reads
        }
        fit = tryCatch(loess(reads ~ covariate, data = x2s, span = 0.3), error = function(e) NULL)
        x$reads = NA
        if (!is.null(fit)) {
          if (is.na(fit$s)) {
            warning("Using all points since initial loess failed")
            fit = loess(reads ~ covariate, data = x2[select, ], span = 1)
          }
        }
        tryCatch(
        {
          if (!is.na(fit$s)) {
            domain = domains[[f]]
            yrange <- quantile(x2s$reads, prob = c(routlier, 1 - routlier), na.rm = TRUE)
            df = data.frame(covariate = seq(domain[1], domain[2], 0.001))
            x$reads = x2$reads/predict(fit, x2) ## apply correction
          }
          else {
            print("Loess failed, yielding NA loess object, continuing without transforming data")
          }
        }, error = function(e) print("Unspecified loess or figure output error"))
      }
      return(x$reads)
    }
  }
  if (verbose) {
    cat('Correcting coverage at individual scales\n')
  }
  ## level 1,2, ..., k corrections
  ## these are the computed corrected values that will be input into the objective function
  correction = NULL
  for (i in rev(1:length(grs))) {
    cat('Correcting coverage at ', WID*base^(i-1), 'bp scale, with', nrow(grs[[i]]), 'intervals\n')
    if (i != length(grs)) {
      grs[[i]]$reads = grs[[i]]$reads/correction[parentmap[[i]][grs[[i]]$lab, parent], cor]
    }
    if (WID*base^(i-1) > 1e5) { ## for very large intervals do not quantile trim, only remove NA


    }
    else {
      grs[[i]]$reads.corrected = FUN(as.data.frame(grs[[i]]), fields, seg = seg);
    }
    if (is.null(correction)) {
      correction = data.table(lab = grs[[i]]$lab, cor = grs[[i]]$reads / grs[[i]]$reads.corrected, key = 'lab')
    }
    else {
      ## multiply new correction and old correction
      old.cor = correction[parentmap[[i]][grs[[i]]$lab, parent], cor]
      new.cor = grs[[i]]$reads / grs[[i]]$reads.corrected                                     
      correction = data.table(lab = grs[[i]]$lab,  cor = old.cor * new.cor, key = 'lab') ## relabel with new nodes
    }
  }
  cov.dt$reads.corrected = grs[[1]][cov.dt$lev0, ]$reads.corrected
  cov.dt[reads.corrected < 0, reads.corrected := NA]
  rm(grs)
  gc()        
  if (verbose) {
    cat('Converting to GRanges\n')
  }
  gc()      
  out = seg2gr(as.data.frame(cov.dt), seqlengths = GenomeInfoDb::seqlengths(cov)) 
  if (verbose) {
    cat('Made GRanges\n')
  }
  gc()
  return(out)
}




#' @name fragCounter
#' @title fragCounter
#' @description Runs entire fragCounter pipeline
#' @author Marcin Imielinski
#' @param bam string path to .bam file
#' @param skeleton string Input data.table with intervals for which there is coverage data
#' @param cov string path to existing coverage rds or bedgraph 
#' @param midpoint boolean If TRUE only count midpoint if FALSE then count bin footprint of every fragment interval
#' @param window integer window / bin size
#' @param minmapq double Minimal map quality
#' @param paired boolean wether or not paired
#' @param outdir string Directory to dump output into
#' @param gc.rds.dir string for tiles of width W, will look here for a file named gc{W}.rds in this directory
#' @param map.rds.dir string for tiles of width W will look here for a file named map{W}.rds in this directory
#' @param exome boolean If TRUE, perform correction using exons as bins instead of fixed size
#' @param use.skel boolean flag If false then default exome skeleton from gencode is used, if TRUE, user defined skeleton is used
#' @export


fragCounter = function(bam, skeleton, cov = NULL, midpoint = TRUE, window = 1000, gc.rds.dir, map.rds.dir, minmapq = 1, paired = TRUE, outdir = NULL, exome = FALSE, use.skel = FALSE) {
  out.rds = paste(outdir, '/cov.rds', sep = '')
  imageroot = gsub('.rds$', '', out.rds)
  if (exome == TRUE) {
    cov = PrepareCov(bam, skeleton = skeleton, cov = NULL, midpoint = midpoint, window = window, minmapq = minmapq, paired = paired, outdir, exome = TRUE, use.skel = use.skel)
    cov = correctcov_stub(cov, gc.rds.dir = gc.rds.dir, map.rds.dir = map.rds.dir, exome = TRUE)
    cov$reads.corrected = coco(cov, mc.cores = 1, fields = c('gc', 'map'), iterative = T, exome = TRUE, imageroot = imageroot)$reads.corrected


  } else {
    cov = PrepareCov(bam, cov = NULL, midpoint = midpoint, window = window, minmapq = minmapq, paired = paired, outdir)
    cov = correctcov_stub(cov, gc.rds.dir = gc.rds.dir, map.rds.dir = map.rds.dir)
    cov$reads.corrected = multicoco(cov, numlevs = 1, base = max(10, 1e5/window), mc.cores = 1, fields = c('gc', 'map'), iterative = T, mono = T)$reads.corrected
  }
  if (!is.null(outdir)) {
    out.rds = paste(outdir, '/cov.rds', sep = '')
    out.corr = paste(gsub('.rds$', '', out.rds), '.corrected.bw', sep = '')
##    if (!is.null(tryCatch({library(rtracklayer); 'success'}, error = function(e) NULL))) { #' twalradt Wednesday, Jan 16, 2019 03:55:28 PM
      cov.corr.out = cov
      cov.corr.out$score = cov$reads.corrected
      cov.corr.out$score[is.na(cov.corr.out$score)] = -1
      cov.corr.out = cov.corr.out[width(cov.corr.out)==window] ## remove any funky widths at end of chromosome
      if (exome == TRUE) {
        export(cov.corr.out[, 'score'], out.corr, 'bigWig', dataFormat = 'variableStep')
      } else {
        export(cov.corr.out[, 'score'], out.corr, 'bigWig', dataFormat = 'fixedStep')
      }
##    } #' twalradt Wednesday, Jan 16, 2019 03:55:58 PM
    saveRDS(cov, paste(gsub('.rds$', '', out.rds), '.rds', sep = ''))
  }
  return(cov)
  cat('done\n')
}




## #' @name bam.cov.exome
## #' @title Get coverage as GRanges from BAM using exons as tiles
## #' @description
## #' Quick way to get tiled coverage via piping to samtools (e.g. ~10 CPU-hours for 100bp tiles, 5e8 read pairs)
## #'
## #' Gets coverage for exons, pulling "chunksize" records at a time and incrementing bin
## #'
## #' @param bam.file string Input BAM file
## #' @param chunksize integer Size of window (default = 1e5)
## #' @param min.mapq integer Minimim map quality reads to consider for counts (default = 30)
## #' @param verbose boolean Flag to increase vebosity (default = TRUE)
## #' @param max.tlen integer Maximum paired-read insert size to consider (default = 1e4)
## #' @param st.flag string Samtools flag to filter reads on (default = '-f 0x02 -F 0x10')
## #' @param fragments boolean flag (default = FALSE) detremining whether to compute fragment (i.e. proper pair footprint aka insert) density or read density
## #' @param do.gc boolean Flag to execute garbage collection via 'gc()' (default = FALSE)
## #' @return GRanges of exon tiles across seqlengths of bam.file with meta data field $counts specifying fragment counts centered (default = TRUE)
## #' in the given bin.
## #' @author Trent Walradt
## #' @export


## bam.cov.exome = function(bam.file, chunksize = 1e6, min.mapq = 1, verbose = TRUE, max.tlen = 1e4, st.flag = "-f 0x02 -F 0x10", fragments = TRUE, do.gc = FALSE)
## {
##   ## check that the BAM is valid
##     check_valid_bam = readChar(gzfile(bam.file, 'r'), 4)
##     if (!identical(check_valid_bam, 'BAM\1')){
##         stop("Cannot open BAM. A valid BAM for 'bam_file' must be provided.")
##     }
##   cmd = 'samtools view %s %s -q %s | cut -f "3,4,9"' ## cmd line to grab the rname, pos, and tlen columns
##   exome = reduce(skidb::read_gencode('exon')) ## Read in exome GRanges to get exons to used as your windows instead of fixed width
##   numwin = length(exome)
##   cat('Calling', sprintf(cmd, st.flag, bam.file, min.mapq), '\n')
##   p = pipe(sprintf(cmd, st.flag, bam.file, min.mapq), open = 'r')
##   i = 0
##   counts = gr2dt(exome)
##   counts[, ':=' (strand = NULL, width = NULL)]
##   setnames(counts, "seqnames", "chr")
##   counts = counts[, bin := 1:length(start)] #, by = chr]
##   counts[, count := 0]
##   counts[, rowid := 1:length(count)]
##   setkeyv(counts, c("chr", "bin")) ## now we can quickly populate the right entries
##   totreads = 0
##   st = Sys.time()
##   if (verbose){
##     cat('Starting fragment count on', bam.file, 'and min mapQ', min.mapq, 'and   insert size limit', max.tlen, '\n')
##   }
##   while (length(chunk <- readLines(p, n = chunksize)) > 0)
##   {
## #  browser(expr = i == 2)
## #  browser()
##     i = i+1
##     chunk = fread(paste(chunk, collapse = "\n"), header = F)[abs(V3) <= max.tlen, ]
##     chunk[, V2 := ifelse(V3 < 0, V2 + V3, V2)] # Convert negative reads to positive
##     chunk[, V3 := abs(V3)]
##     if (grepl("chr",chunk$V1[1])) { # Robust to samples that start with or without 'chr' before chromosomes
##       chunk[, V1 := gsub("chr", "", V1)]
##     }
##     chunk = chunk[which(V1 %in% GenomeInfoDb::seqlevels(exome))] ## Exome only has seqlevels 1-22,X,Y,M, remove any additional seqlevels from sample
##     if (nrow(chunk) > 0) {
##       chunk.gr = GRanges(seqnames = chunk$V1, ranges = IRanges(start = chunk$V2, width = chunk$V3))
##       ## Robust to chunks that fall entirely between exons
##       chunk.match = tryCatch(
##         gr.match(chunk.gr,exome),
##         error = function(e) e
##       )
##       if(!inherits(chunk.match, "error")){
##         chunk[, bin := chunk.match]
##         tabs = chunk[, list(newcount = length(V1)), by = list(chr = as.character(V1), bin)] ## tabulate reads to bins data.table style
##         counts[tabs, count := count + newcount] ## populate latest bins in master data.table
##       }
##     }
##     ## should be no memory issues here since we preallocate the data table .. but they still appear
##     if (do.gc){
##         print('GC!!')
##         print(gc())
##     }
##     ## report timing
##     if (verbose){
##       cat('bam.cov.tile.st ', bam.file, 'chunk', i, 'num fragments processed', i*chunksize, '\n')
##       timeelapsed = as.numeric(difftime(Sys.time(), st, units = 'hours'))
##       meancov = i * chunksize / counts[tabs[nrow(tabs),], ]$rowid  ## estimate comes from total reads and "latest" bin filled
##       totreads = meancov * numwin
##       tottime = totreads*timeelapsed/(i*chunksize)
##       rate = i*chunksize / timeelapsed / 3600
##       cat('mean cov:', round(meancov,1), 'per bin, estimated tot fragments:', round(totreads/1e6,2), 'million fragments, processing', rate,
##           'fragments/second\ntime elapsed:', round(timeelapsed,2), 'hours, estimated time remaining:', round(tottime - timeelapsed,2), 'hours', ', estimated total time', round(tottime,2), 'hours\n')
##     }
##   }
##   x = data.table(chr = c(1:22, "X", "Y", "M"), order = 1:25)
##   counts[, order := x$order[match(chr, x$chr)]]
##   counts = counts[order(order)][, order := NULL]
##   exome$counts = counts$count
##     if (verbose){
##       cat("Finished computing coverage, and making GRanges\n")
##     }
##   close(p)
##   return(exome)
## }




#' @title Mappability calculator
#' @description Calculates mappability as fraction of bases in a tile with mappability of 1
#' @name MAP.fun
#' @param win.size integer Window size, in basepairs, to calculate mappability for (should match window
#' size of your coverage file
#' @param twobitURL string URL to twobit genome file. Default is hg19 from UCSC
#' @param bw.path string Path to .bigWig mappability file
#' @param twobit.win integer How many windows of the twobit file to load into memory on each core
#' @param mc.cores integer How many cores to use
#' @param exome boolean If TRUE, calculate mappability for exons instead of window
#' @author Trent Walradt
#' @export


MAP.fun = function(win.size = 200, twobitURL = '~/DB/UCSC/hg19.2bit', bw.path = '~/DB/UCSC/wgEncodeCrgMapabilityAlign100mer.bigWig', twobit.win = 1e3, mc.cores = 1, exome = FALSE) {
  if (exome == TRUE){
    tiles = reduce(read_gencode('exon'))
    GenomeInfoDb::seqlevelsStyle(tiles) <- "UCSC" #Formatting chromosome notation
  } else {
    tiles = gr.tile(GenomeInfoDb::seqlengths(TwoBitFile(twobitURL)), win.size)
  }
  values(tiles) = NULL
  tiles = tiles %Q% (seqnames %in% c(paste0("chr",seq(1,22)),"chrX","chrY")) # removes all chromomes except 1:22 and X/Y
  GenomeInfoDb::seqlevels(tiles) = GenomeInfoDb::seqlevelsInUse(tiles) # Get rid of empty seqname factor levels
  x = seq(1,(length(tiles) / twobit.win)) 
  map.score = function(x) {
    this.ix = seq(1 + ((x-1) * twobit.win), x * twobit.win)
    out = data.table(mappability = tiles[this.ix] %O% (rtracklayer::import(bw.path, selection = tiles[this.ix]) %Q% (score==1)), index = this.ix)
    return(out)
  }
  map.out = mclapply(x, map.score, mc.cores = mc.cores)
  map.out = rbindlist(map.out)
  if (!is.integer(length(tiles)/twobit.win)){
    edge.num = seq(twobit.win*max(x)+1,length(tiles))
    edge.out = data.table(mappability = tiles[edge.num] %O% (rtracklayer::import(bw.path, selection = tiles[edge.num]) %Q% (score==1)), index = edge.num)
    map.out = rbind(map.out, edge.out)
  }
  setkey(map.out,index)
  tiles$score = map.out$mappability
  tiles = gr.sub(tiles)
  tiles = gr.stripstrand(tiles)
  return(tiles)
}




#' @title GC content calculator
#' @description Calculates GC content across the genome for a given sized window
#' @name GC.fun
#' @param win.size integer Window size, in basepairs, to calculate GC content for (should match window
#' size of your coverage file
#' @param twobitURL string URL to twobit genome file. Default is hg19 from UCSC
#' @param twobit.win integer How many windows of the twobit file to load into memory on each core
#' @param mc.cores integer How many cores to use
#' @param exome boolean If TRUE, calculate mappability for exons instead of window
#' @author Trent Walradt
#' @export


GC.fun = function(win.size = 200, twobitURL = '~/DB/UCSC/hg19.2bit', twobit.win = 1e3, mc.cores = 1, exome = FALSE) {
  if (exome == TRUE){
    tiles = reduce(read_gencode('exon'))
    GenomeInfoDb::seqlevelsStyle(tiles) <- "UCSC"
  } else {
    tiles = gr.tile(GenomeInfoDb::seqlengths(TwoBitFile(twobitURL)),win.size)
  }
  values(tiles) = NULL
  tiles = tiles %Q% (seqnames %in% c(paste0("chr",seq(1,22)),"chrX","chrY")) # removes all chromomes except 1:22 and X/Y
  x = seq(1,(length(tiles) / twobit.win))
  gc.con = function(x) {
    this.ix = seq(1 + ((x-1) * twobit.win), x * twobit.win)
    # print(this.ix[1])
    tmp = alphabetFrequency(ffTrack::get_seq(twobitURL, tiles[this.ix]))
    out = data.table(gc = rowSums(tmp[, c('C', 'G')])/rowSums(tmp), index = this.ix)
    return(out)
  }
  gc.out = mclapply(x, gc.con, mc.cores = mc.cores)
    gc.out = rbindlist(gc.out)
  if (!is.integer(length(tiles)/twobit.win)){
    edge.num = seq(twobit.win*max(x)+1,length(tiles))
    tmp.edge = alphabetFrequency(ffTrack::get_seq(twobitURL, tiles[edge.num]))
    edge.out = data.table(gc = rowSums(tmp.edge[, c('C', 'G')])/rowSums(tmp.edge), index = edge.num)
    gc.out = rbind(gc.out, edge.out)
  }
  setkey(gc.out,index)
  tiles$score = gc.out$gc
  tiles = gr.sub(tiles)
  tiles = gr.stripstrand(tiles)
  return(tiles)
}




#' @name PrepareCov
#' @title PrepareCov
#' @description Load BAM or coverage file and prepare for use in correctcov_stub
#' @author Marcin Imielinski
#' @param bam path to .bam file
#' @param skeleton string Input data.table with intervals for which there is coverage data
#' @param cov Path to existing coverage rds or bedgraph 
#' @param midpoint If TRUE only count midpoint if FALSE then count bin footprint of every fragment interval
#' @param window window / bin size
#' @param minmapq Minimal map quality
#' @param paired wether or not paired
#' @param outdir Directory to dump output into
#' @param exome boolean If TRUE, use bam.cov.exome to calculate coverage
#' @param use.skel boolean flag If false then default exome skeleton from gencode is used, if TRUE, user defined skeleton is used
#' @author Trent Walradt
#' @export


PrepareCov = function(bam, skeleton, cov = NULL, midpoint = TRUE, window = 200, minmapq = 1, paired = TRUE, outdir = NULL, exome = FALSE, use.skel = FALSE) {
  if (exome == TRUE){
#    cov = bam.cov.exome(bam, chunksize = 1e6, min.mapq = 1)
    cov = bam.cov.skel(bam, skeleton, chunksize = 1e6, min.mapq = 1, use.skel = use.skel)
  } else {
    if (is.null(bam)) {
      bam = ''
    }


    if (file.exists(bam)) { # & is.null(cov))
      if (!midpoint) {
        cat("Running without midpoint!!!\n")
      }
      print('Doing it!')
      if (is.null(paired)) {
        paired = TRUE
      }
      if (paired) {
        cov = bamUtils::bam.cov.tile(bam, window = as.numeric(window), chunksize = 1e6, midpoint = TRUE, min.mapq = 1)  ## counts midpoints of fragments
      }
      else {
        sl = GenomeInfoDb::seqlengths(BamFile(bam))
        tiles = gr.tile(sl, window)
        cov = bamUtils::bam.cov.gr(bam, intervals = tiles, isPaired = NA, isProperPair = NA, hasUnmappedMate = NA, chunksize = 1e5, verbose = TRUE)  ## counts midpoints of fragments    # Can we increase chunksize?
        cov$count = cov$records/width(cov)
      }
    }
    else if (!is.null(cov)) {
      cov = readRDS(cov)
    }
    else {
      stop("Can't locate either bam or coverage input file")
    }
  }
  gc()
  cat('Finished acquiring coverage\n')
  return(cov)
  cat('done\n')
}




#' @title Correct coverage stub
#' @description prepares GC, mappability, and coverage files for multicoco
#' @name correctcov_stub
#' @param cov.wig wig file of coverage tiles of width W or pointer to rds file
#' of sorted GRanges object or GRanges object
#' @param mappability double threshold for mappability score
#' @param samplesize integer size of sub-sample
#' @param verbose boolean Wether to print log to console
#' @param gc.rds.dir string for tiles of width W, will look here for a file named gc{W}.rds in this directory
#' @param map.rds.dir string for tiles of width W will look here for a file named map{W}.rds in this directory
#' @param exome boolean If TRUE, look in gc/map.rds.dir for files called 'gcexome.rds' and 'mapexome.rds'
#' @author Trent Walradt
#' @export


correctcov_stub = function(cov.wig, mappability = 0.9, samplesize = 5e4, verbose = T, gc.rds.dir, map.rds.dir, exome = FALSE) {
  if (is.character(cov.wig)) {
    if (grepl('(\\.bedgraph$)|(\\.wig$)|(\\.bw$)', cov.wig)) {
      cov = import.ucsc(cov.wig)
    }
    else if (grepl('(\\.rds)', cov.wig)) {
      cov = gr.stripstrand(readRDS(cov.wig))
      names(values(cov))[1] = 'score' ## will take the first value column as the (uncorrected) read count
    }
    else {
      stop("Unsupported coverage format (should be UCSC bedgraph, wig, bw, or R Bioconductor GRanges .rds file")
    }
  }
  else { ## assume it is a sorted GRanges
    cov = gr.stripstrand(cov.wig)
  }
  # cov = cov %Q% (seqnames == 21) #' twalradt Monday, Apr 23, 2018 10:33:19 AM Done for unit testing
  n = length(cov)
  wid = as.numeric(names(sort(-table(width(cov))))[1])
  if (exome == TRUE) {
    gc.rds = paste(gc.rds.dir,'/gcexome.rds', sep = '')
    map.rds = paste(map.rds.dir,'/mapexome.rds', sep = '')
  } else {
    gc.rds = paste(gc.rds.dir,'/gc', wid, '.rds', sep = '')
    map.rds = paste(map.rds.dir,'/map', wid, '.rds', sep = '')
  }
  cat('Loaded GC and mappability\n')
  if (!file.exists(gc.rds) | !file.exists(map.rds)) {
    stop(sprintf('GC rds file %s not found, either directory is wrong or file needs to be generated for this particular window width', gc.wig))
    stop(sprintf('mappability rds file %s not found, either directory is wrong or file needs to be generated for this particular window width', map.wig))
  }
  else { ## if we have rds files then let's use these to avoid using rtracklayer
    gc = readRDS(gc.rds)
    map = readRDS(map.rds)
  }
  if (is.null(cov$score)) { ## if $score field is empty then just assume that the first column of coverage is the "score" i.e. read count
    names(values(cov))[1] = 'score'
  }
  #map = gr.sub(map)
  #gc = gr.sub(gc)
  gc.str = gr.string(gc)
  map.str = gr.string(map)
  cov.str = gr.string(cov)
  all.str = intersect(intersect(map.str, cov.str), gc.str) ## in case we have mismatches in the ordering / genome definition       
  cat(sprintf('length cov is %s, length gc is %s, length map is %s\n',
              length(cov),
              length(gc),
              length(gc)))
  map = map[match(all.str, map.str)]
  cov = cov[match(all.str, cov.str)]
  gc = gc[match(all.str, gc.str)]        
  if (length(cov) != length(gc) | length(gc) != length(map)) {
    stop('Mismatch / problem in cov, gc, or map definition.  Check if they come from the same width tiling')
  }  
  cov$reads = cov$score
  cov$gc = gc$score
  cov$gc[cov$gc<0] = NA
  cov$map = map$score
  cov$score = NULL
  rm(gc)
  rm(map)
  gc()  
  cat('Synced coverage, GC, and mappability\n')
  cov = sort(gr.fix(cov))        
  cat('Modified gc / mappability correction\n')
  return(cov)
}






#' @title Coverage correction
#' @description Given gc and mappability coverage correction, yields corrected
#' read counts
#' @name coco
#' @param cov GRanges constant with GRanges of coverage samples with (by default) fields $reads, $map, $gc
#' @param base integer Scale multiplier
#' @param fields character vector fields of gc to use as covariates
#' @param iterative boolean whether to iterate
#' @param presegment boolean whether to presegment
#' @param min.segwidth integer when presegmenting, minimum segment width
#' @param verbose boolean Wether to print log to console
#' @param FUN function with which to correct coverage (by default loess
#' correction modified from HMMcopy that takes in granges with fields
#' $reads and other fields specified in "fields"
#' @param ... additional args to FUN
#' @param mc.cores integer Number of cores to use
#' @param exome boolean If TRUE, collapse by 1e5 for presegmentation
#' @param imageroot String, optional file root to which to dump images of correction
#' @author Trent Walradt
#' @usage coco(cov, base = max(10, 1e5/max(width(cov))), fields = c("gc", "map"), iterative = TRUE,
#' presegment = TRUE, min.segwidth = 5e6, verbose = TRUE,
#' FUN = NULL, ..., mc.cores = 1, exome = TRUE, imageroot = NULL)
#' @export


coco = function(cov, base = max(10, 1e5 / max(width(cov))), fields = c("gc", "map"),
                     iterative = TRUE, presegment = TRUE, min.segwidth = 5e6, verbose = TRUE, FUN = NULL, ..., mc.cores = 1, exome = TRUE, imageroot = NULL) {
  if (verbose) {
    cat('Converting to data.table\n')
  }


  cov = sort(sortSeqlevels(cov))
  WID = max(width(cov))
##  library(data.table) #' twalradt Wednesday, Jan 16, 2019 03:54:02 PM 
  cov.dt = gr2dt(cov)        
  sl = structure(as.numeric(1:length(GenomeInfoDb::seqlevels(cov))), names = GenomeInfoDb::seqlevels(cov))       
  if (verbose) {
    cat('Grouping intervals\n')
  }
  ##  -.01 just to make Trent Happy ie produce identical result to above
  cov.dt[, lev1 :=  as.character(sl[seqnames] + ceiling((1:.N-0.01)/base)), by = seqnames]
  if (presegment) { ## perform rough segmentation at highest level
    seg = NULL
    sl = GenomeInfoDb::seqlengths(cov)
    if (verbose) {
      cat('Presegmenting at ', as.integer(WID*base), ' bp scale \n')
    }
##    require(DNAcopy) #' twalradt Wednesday, Jan 16, 2019 04:05:49 PM
    set.seed(42) #' twalradt Friday, Apr 20, 2018 01:07:28 PM
    if (exome == TRUE) {
      tiles = gr.tile(sl, width = 1e5)
      tiles = tiles %Q% (seqnames %in% c(seq(1,22),"X","Y"))
      GenomeInfoDb::seqlevels(tiles) = GenomeInfoDb::seqlevelsInUse(tiles) # Get rid of empty seqname factor levels
      tmp.cov = tiles %$% cov
      tmp.cov = tmp.cov %Q% (!is.na(reads))
      # tmp.cov = gr.val(tiles, cov, val = c('reads', 'gc', 'map'), na.rm = TRUE)
    } else {
      tmp.cov = seg2gr(cov.dt[,list(chr = seqnames[1], start = min(start), end = max(end), strand = strand[1], reads = mean(reads, na.rm = T)), by = get(paste("lev", numlevs, sep = ''))][end>start, ], seqlengths = sl)
    }
    ix = which(!is.na(values(tmp.cov)[, 'reads']))
    tmp = data.frame()
    tryCatch({
      cna = CNA(log(values(tmp.cov)[, 'reads'])[ix], as.character(seqnames(tmp.cov))[ix], start(tmp.cov)[ix], data.type = 'logratio')
      tmp = print(segment(smooth.CNA(cna), alpha = 1e-5, verbose = T))
      tmp = tmp[!is.na(tmp$loc.start) & !is.na(tmp$chrom) & !is.na(tmp$loc.end), , drop = F]
    }, error = function(e) warning('DNACopy error moving on without segmenting'))
    if (nrow(tmp)>0) {
      seg = sort(seg2gr(tmp, seqlengths = sl))
      seg = seg[width(seg)>min.segwidth] ## remove small segments
      seg.dt = gr2dt(seg)
      if (nrow(seg.dt)>0) {
        seg = seg2gr(seg.dt[, list(seqnames = seqnames,
                                   start = ifelse(c(FALSE, seqnames[-length(seqnames)]==seqnames[-1]), c(1, start[-1]), 1),
                                   end = ifelse(c(seqnames[-length(seqnames)]==seqnames[-1], FALSE), c(start[-1]-1, Inf), GenomeInfoDb::seqlengths(seg)[as.character(seqnames)]))], seqlengths = sl)
        seg = gr.val(seg, tmp.cov, 'reads') ## populate with mean coverage
        seg$reads = seg$reads/sum(as.numeric(seg$reads*width(seg))/sum(as.numeric(width(seg)))) ## normalize segs by weigthed mean (so these become a correction factor)
      }
      else {
        seg = NULL
      }
    }
    else {
      seg = NULL
    }
  }
  else {
    seg = NULL
  }
  ## modified from HMMCopy to
  ## (1) take arbitrary set of covariates, specified by fields vector
  ## (2) employ as input an optional preliminary (coarse) segmentation with which to adjust signal immediately prior to loess
  ## NOTE: (this only impacts the loess fitting, does not impose any segmentation on the data)
  ##
  if (is.null(FUN)) {
    FUN = function(x, fields = fields, samplesize = 5e4, seg = NULL, verbose = T, doutlier = 0.001, routlier = 0.01) {
      if (!all(fields %in% names(x))) {
        stop(paste('Missing columns:', paste(fields[!(fields %in% names(x))], collapse = ',')))
      }
      x$valid <- TRUE
      x = as.data.frame(x)
      for (f in fields) {
        x$valid[is.na(x[, f])] = FALSE
        x$valid[which(is.infinite(x[, f]))] = FALSE
      }
      if (verbose) {
        cat('Quantile filtering response and covariates\n')
      }
      range <- quantile(x$reads[x$valid], prob = c(routlier, 1 - routlier), na.rm = TRUE)
      if (verbose) {
        cat(sprintf("Response min quantile: %s max quantile: %s\n", round(range[1],2), round(range[2],2)))
      }
      domains = lapply(fields, function(f) quantile(x[x$valid, f], prob = c(doutlier, 1 - doutlier), na.rm = TRUE))
      names(domains) = fields
      x$ideal <- x$valid
      x$ideal[x$reads<=range[1] | x$reads>range[2]] = FALSE
      for (f in fields) {
        x$ideal[x[, f] < domains[[f]][1] | x[, f] > domains[[f]][2]] = FALSE
      }
      if (verbose) {
        cat(sprintf('Nominated %s of %s data points for loess fitting\n', sum(x$ideal), nrow(x)))
      }
      set <- which(x$ideal)
      if (length(set)<=10) {
        warning("Not enough samples for loess fitting - check to see if missing or truncated data?")
        return(x$reads)
      }
      for (f in fields) {
        if (verbose) {
          message(sprintf("Correcting for %s bias...", f))
        }
        set.seed(42)
        select <- sample(set, min(length(set), samplesize))
        x2 = x[, c('reads', f)]
        x2$covariate = x2[, f]
        x2s = x2[select, ]                    
        if (!is.null(seg)) {  ## here we apply a prelmiinary segmentation to correct for large scale copy variation allow more power to reveal the covariate signal
          if (verbose) {  
            message('Applying preliminary segmentation prior to loess fitting')
          }
          x.grs = gr.val(seg2gr(x[select, ], seqlengths = NULL), seg, 'reads')
          x2s$reads = x2s$reads/x.grs$reads 
        }
        fit = tryCatch(loess(reads ~ covariate, data = x2s, span = 0.3), error = function(e) NULL)
        x$reads = NA
        if (!is.null(fit)) {
          if (is.na(fit$s)) {
            warning("Using all points since initial loess failed")
            fit = loess(reads ~ covariate, data = x2[select, ], span = 1)
            domain = domains[[f]]
            yrange <- quantile(x2s$reads, prob = c(routlier, 1 - routlier), na.rm = TRUE)
            df = data.frame(covariate = seq(domain[1], domain[2], 0.001))
            
            if (!is.null(imageroot)) {
              out.png = paste(imageroot, ifelse(grepl("/$", imageroot), '', '.'), f,'_correction.png', sep = '')
              if (verbose) {
                cat("Dumping figure to", out.png, "\n")
              }
              png(out.png, height = 1000, width = 1000) 
              plot(x2s$covariate, x2s$reads, col = alpha('black', 0.1), pch = 19, cex = 0.4, xlim = domain, ylim = yrange, ylab = sprintf('signal before %s correction', f), xlab = f);
              lines(df$covariate, predict(fit, df), col = 'red', lwd = 2)
              dev.off()
            }
          }
        }
        tryCatch(
        {
          if (!is.na(fit$s)) {
            domain = domains[[f]]
            yrange <- quantile(x2s$reads, prob = c(routlier, 1 - routlier), na.rm = TRUE)
            df = data.frame(covariate = seq(domain[1], domain[2], 0.001))


            if (!is.null(imageroot)) {
              out.png = paste(imageroot, ifelse(grepl("/$", imageroot), '', '.'), f,'_correction.png', sep = '')
              if (verbose) {
                cat("Dumping figure to", out.png, "\n")
              }
              png(out.png, height = 1000, width = 1000) 
              plot(x2s$covariate, x2s$reads, col = alpha('black', 0.1), pch = 19, cex = 0.4, xlim = domain, ylim = yrange, ylab = sprintf('signal before %s correction', f), xlab = f);
              lines(df$covariate, predict(fit, df), col = 'red', lwd = 2)
              dev.off()
            }


            
            predicted.fit = predict(fit, x2) ## apply correction
            reads.x2 = x2$reads
            x$reads = x2$reads/predicted.fit ## Split this into 2 lines


          }
          else {
            print("Loess failed, yielding NA loess object, continuing without transforming data")
          }
        }, error = function(e) print("Unspecified loess or figure output error"))
      }
      return(list(x$reads, predicted.fit, reads.x2))
    }
  }
  FUN.out = FUN(as.data.frame(cov.dt), fields, seg = seg)
  cov.dt$reads.corrected = FUN.out[[1]]
#browser()
  cov.dt[reads.corrected < 0, reads.corrected := NA] #' twalradt Tuesday, Oct 09, 2018 11:46:47 AM  ADDED THIS LINE and commented out the next line
  if (verbose) {
    cat('Converted negatives to Nas\n')
  }
                                        #  cov.dt[, reads.corrected := pmax(0,reads.corrected)]
  cov.dt$predicted.fit = FUN.out[[2]]
  cov.dt$reads.x2 = FUN.out[[3]]


  #cov.dt$reads.corrected = FUN(as.data.frame(cov.dt), fields, seg = seg)[[1]] # Returning first item (x$reads)
  gc()        
  if (verbose) {
    cat('Converting to GRanges\n')
  }
  gc()      
  out = seg2gr(as.data.frame(cov.dt), seqlengths = GenomeInfoDb::seqlengths(cov)) 
  if (verbose) {
    cat('Made GRanges\n')
  }
  gc()
  return(out)
}








#' @title alpha
#' @description
#' Takes provided colors and gives them the specified alpha (ie transparency) value
#' @name alpha
#' @param col RGB color
#' @param alpha value of alpha


alpha = function(col, alpha)
{    
  col.rgb = col2rgb(col)
  out = rgb(red = col.rgb['red', ]/255, green = col.rgb['green', ]/255, blue = col.rgb['blue', ]/255, alpha = alpha)
  names(out) = names(col)
  return(out)
}








#' @name bam.cov.skel
#' @title Get coverage as GRanges from BAM using exons as tiles
#' @description
#' Quick way to get tiled coverage via piping to samtools (e.g. ~10 CPU-hours for 100bp tiles, 5e8 read pairs)
#'
#' Gets coverage for user-defined regions of the genome, pulling "chunksize" records at a time and incrementing bin
#'
#' @param bam.file string Input BAM file
#' @param skeleton string Input data.table with intervals for which there is coverage data
#' @param chunksize integer Size of window (default = 1e5)
#' @param min.mapq integer Minimim map quality reads to consider for counts (default = 30)
#' @param verbose boolean Flag to increase vebosity (default = TRUE)
#' @param max.tlen integer Maximum paired-read insert size to consider (default = 1e4)
#' @param st.flag string Samtools flag to filter reads on (default = '-f 0x02 -F 0x10')
#' @param fragments boolean flag (default = FALSE) detremining whether to compute fragment (i.e. proper pair footprint aka insert) density or read density
#' @param do.gc boolean Flag to execute garbage collection via 'gc()' (default = FALSE)
#' @param use.skel boolean flag If false then default exome skeleton from gencode is used, if TRUE, user defined skeleton is used
#' @return GRanges of user defined tiles across seqlengths of bam.file with meta data field $counts specifying fragment counts centered (default = TRUE)
#' in the given bin.
#' @author Trent Walradt
#' @export


bam.cov.skel = function(bam.file, skeleton, chunksize = 1e5, min.mapq = 1, verbose = TRUE, max.tlen = 1e4, st.flag = "-f 0x02 -F 0x10", fragments = TRUE, do.gc = FALSE, use.skel = FALSE)
{
  ## check that the BAM is valid
    check_valid_bam = readChar(gzfile(bam.file, 'r'), 4)
    if (!identical(check_valid_bam, 'BAM\1')){
        stop("Cannot open BAM. A valid BAM for 'bam_file' must be provided.")
    }
  cmd = 'samtools view %s %s -q %s | cut -f "3,4,9"' ## cmd line to grab the rname, pos, and tlen columns
  
  # numwin = length(exome)


  if (use.skel == FALSE){
    skel = reduce(skidb::read_gencode('exon'))
    counts = gr2dt(skel)
  } else {
    skel = readRDS(skeleton)
    counts = skel
  }


  numwin = nrow(skel)
  cat('Calling', sprintf(cmd, st.flag, bam.file, min.mapq), '\n')
  p = pipe(sprintf(cmd, st.flag, bam.file, min.mapq), open = 'r')
  i = 0


  # counts = gr2dt(exome)
  counts[, ':=' (strand = NULL, width = NULL)]
  setnames(counts, "seqnames", "chr")
  counts = counts[, bin := 1:length(start)] #, by = chr]
  counts[, count := 0]
  counts[, rowid := 1:length(count)]
  counts$chr <- factor(counts$chr, levels = c(as.character(1:22), "X", "Y")) #Ensure counts is sorted with chromsomes going 1:22,x,y
  setkeyv(counts, c("chr", "bin")) ## now we can quickly populate the right entries
  totreads = 0
  st = Sys.time()
  if (verbose){
    cat('Starting fragment count on', bam.file, 'and min mapQ', min.mapq, 'and   insert size limit', max.tlen, '\n')
  }
  while (length(chunk <- readLines(p, n = chunksize)) > 0)
  {
#  browser(expr = i == 86)
    i = i+1
    chunk = fread(paste(chunk, collapse = "\n"), header = F)[abs(V3) <= max.tlen, ]
    chunk[, V2 := ifelse(V3 < 0, V2 + V3, V2)] # Convert negative reads to positive
    chunk[, V3 := abs(V3)]
    if (grepl("chr",chunk$V1[1])) { # Robust to samples that start with or without 'chr' before chromosomes
      chunk[, V1 := gsub("chr", "", V1)]
    }
    if (use.skel == FALSE){ chunk = chunk[which(V1 %in% GenomeInfoDb::seqlevels(skel))]
    } else { chunk = chunk[which(V1 %in% levels(skel$chr))]
    }
    if (nrow(chunk) > 0) {
      chunk.gr = GRanges(seqnames = chunk$V1, ranges = IRanges(start = chunk$V2, width = chunk$V3))
      ## Robust to chunks that fall entirely between exons
      if (use.skel == TRUE) {
        chunk.match = tryCatch(
          gr.match(chunk.gr,dt2gr(skel)),
          error = function(e) e
        )
      } else {
        chunk.match = tryCatch(
          gr.match(chunk.gr,skel),
          error = function(e) e
        )
      }


      if(!inherits(chunk.match, "error")){
        chunk[, bin := chunk.match]
        tabs = chunk[, list(newcount = length(V1)), by = list(chr = as.character(V1), bin)] ## tabulate reads to bins data.table style
        counts[tabs, count := count + newcount] ## populate latest bins in master data.table
      }
    }
    ## should be no memory issues here since we preallocate the data table .. but they still appear
    if (do.gc){
        print('GC!!')
        print(gc())
    }
    ## report timing
    if (verbose){
      cat('bam.cov.tile.st ', bam.file, 'chunk', i, 'num fragments processed', i*chunksize, '\n')
      timeelapsed = as.numeric(difftime(Sys.time(), st, units = 'hours'))
      meancov = i * chunksize / counts[tabs[nrow(tabs),], ]$rowid  ## estimate comes from total reads and "latest" bin filled
      totreads = meancov * numwin
      tottime = totreads*timeelapsed/(i*chunksize)
      rate = i*chunksize / timeelapsed / 3600
      cat('mean cov:', round(meancov,1), 'per bin, estimated tot fragments:', round(totreads/1e6,2), 'million fragments, processing', rate,
          'fragments/second\ntime elapsed:', round(timeelapsed,2), 'hours, estimated time remaining:', round(tottime - timeelapsed,2), 'hours', ', estimated total time', round(tottime,2), 'hours\n')
    }
  }
  x = data.table(chr = c(1:22, "X", "Y", "M"), order = 1:25)
  counts[, order := x$order[match(chr, x$chr)]]
  counts = counts[order(order)][, order := NULL]


  if (use.skel == TRUE){
    skel$counts = counts$count
    skel = dt2gr(skel)
  } else {
    skel$counts = counts$count
  }  
  if (verbose){
    cat("Finished computing coverage, and making GRanges\n")
  }
  close(p)
  return(skel)
}




#' @name make.blacklist
#' @title Make a list of regions of the genome in which there is coverage in the normal samples
#' @description
#' Look at all available normal samples and find regions where some percentage of all normal samples have at least 1 read
#'
#' @param pairs string Input data.table with columns 'pair' and 'bam' corresponding to the sample name and bath to normal bam
#' @param cutoff double Threshold for what percentage of normal samples must have reads for the region to be included
#' @param mc.cores integer Number of cores to use
#' @return Data.table of regions to include when calculating corrected reads
#' @author Trent Walradt
#' @export


make.blacklist = function(pairs, cutoff = 0.9, mc.cores = mc.cores)
{


  pairs.bl = readRDS(pairs)
  setkeyv(pairs.bl, "pair")


  bl = mclapply(pairs.bl[, pair], function(nm){
    this.cov = bam.cov.skel(pairs.bl[nm, bam], chunksize = 1e6)
    if (!is.null(this.cov)){
      this.cov = transpose(gr2dt(tmp)[,.(counts)])
      message(nm)
      print(dim(this.cov))
    } else {this.cov = data.table(NA)}
    return(this.cov)}
  , mc.cores = mc.cores)


  bl.bind = rbindlist(bl, fill = T)
  bl.bind.t = transpose(bl.bind)
  for(col in names(bl.bind.t)) set(bl.bind.t, i = which(bl.bind.t[[col]] < 1 ), j = col, value = NA)
  for(col in names(bl.bind.t)) set(bl.bind.t, i = which(!is.na(bl.bind.t[[col]])), j = col, value = 1)
  for(col in names(bl.bind.t)) set(bl.bind.t, i = which(is.na(bl.bind.t[[col]])), j = col, value = 0)
  bl.bind.t[, black_list_pct := rowSums(.SD)/dim(bl.bind.t)[2]]
  bl.bind.t[, blacklisted := ifelse(black_list_pct > cutoff, "Keep", "Remove")]
  exome = gr2dt(reduce(read_gencode('exon')))[, blacklisted := bl.bind.t$blacklisted][blacklisted == "Keep"][, blacklisted := NULL]
  return(exome)


}













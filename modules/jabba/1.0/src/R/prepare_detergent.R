#' @import GenomicRanges
#' @importFrom data.table data.table
#' @importFrom data.table rbindlist
#' @importFrom data.table set
#' @importFrom data.table setDT
#' @importFrom data.table setkey
#' @importFrom data.table setkeyv
#' @importFrom data.table setnames
#' @importFrom data.table transpose
#' @import rsvd
#' @importFrom magrittr %>%
#' @importFrom gUtils gr2dt
#' @importFrom gUtils dt2gr
#' @importFrom stats hclust
#' @importFrom stats cutree
#' @importFrom stats dist
#' @importFrom stats median
#' @importFrom stats na.omit
#' @importFrom MASS ginv
#' @importFrom utils globalVariables
#' @import gUtils
#' @import pbmcapply


suppressPackageStartupMessages({
    library(GenomicRanges)
    library(data.table)
    library(rsvd)
    library(magrittr)
    library(gUtils)
    library(stats)
    library(MASS)
    library(utils)
    library(pbmcapply)
})

globalVariables(c(".", "..ix", "L", "L1", "V1", "black_list_pct", "blacklisted", "decomposed_cov", "germline.status", "log.reads", "mclapply", "median.chr", "normal_cov", "foreground", "input.read.counts", "foreground.log", "reads.corrected", "background", "file.available", "mt", "background.log", ".N", ".SD", ":=", "median.idx", ".GRP", "reads.corrected.org", "%>%", "signal", "signal.org"))

prepare_detergent <- function(normal.table.path = NA, use.all = TRUE, choose.randomly = FALSE, choose.by.clustering = FALSE, number.of.samples = 50, save.pon = FALSE, path.to.save = NA, verbose = TRUE, num.cores = 1, tolerance = 0.0001, is.human = TRUE, build = "hg19", field = "reads.corrected", PAR.file = NULL){
    
    if (verbose){
        message("Starting the preparation of Panel of Normal samples a.k.a detergent")
    }


    if (is.na(normal.table.path)){
        stop("Need a table with paths to normal samples to create a PON")
    }


    if (is.na(path.to.save) & save.pon == TRUE){
        stop("Need a path to save decomposed PON")
    }


    normal.table = readRDS(normal.table.path)
    setkeyv(normal.table, "sample")


    num.samp = nrow(normal.table)


    if (verbose){
        message(paste0(num.samp, " samples available"))
    }
        
    if (use.all & choose.randomly | use.all & choose.by.clustering | choose.randomly & choose.by.clustering | use.all & choose.randomly & choose.by.clustering){
        stop("only one of use.all, choose.randomly, choose.by.clustering can be set to TRUE. Rectify and restart")
    }


    if (choose.randomly){
        if (verbose){ 
            message(paste0("Selecting ", number.of.samples, " normal samples randomly"))
        }
        set.seed(12)
        samp.final = sample(1:num.samp, number.of.samples)
        samp.final = normal.table[samp.final]
        setkey(samp.final, "sample")
    }


    if (choose.by.clustering){


        if (verbose){ 
            message("Starting the clustering")
        }
        
        mat.small = mclapply(normal.table[, sample], function(nm){
            this.cov = tryCatch(readRDS(normal.table[nm, normal_cov]), error = function(e) NULL)
            if (!is.null(this.cov)){
                this.cov = this.cov[, field] %>% gr2dt() %>% setnames(., field, "signal")
                ## reads = this.cov[seqnames == "22", .(seqnames, signal)]
                reads = this.cov[seqnames == seqnames[1], .(seqnames, signal)]
                reads[, median.chr := median(.SD$signal, na.rm = T), by = seqnames]
                reads[is.na(signal), signal := median.chr]
                min.cov = min(reads[signal > 0]$signal, na.rm = T)
                reads[signal == 0, signal := min.cov]
                reads[signal < 0, signal := min.cov]
                reads = log(reads[, .(signal)])
                reads = transpose(reads)
                reads = cbind(reads, nm)
            } else {reads = data.table(NA)}
            return(reads)}, mc.cores = num.cores)
        
        gc()
        
        mat.sub = rbindlist(mat.small, fill = T)
        mat.sub = na.omit(mat.sub)
        ix = ncol(mat.sub)
        samp.names = mat.sub[, ..ix]
        mat.sub.t = transpose(mat.sub[, 1:(ncol(mat.sub) - 1)])
        rm(mat.sub)
        gc()


        if (verbose){ 
            message("Starting decomposition on a small section of genome")
        }


        mat.sub.t = as.matrix(mat.sub.t)
        gc()


        rpca.mat = rrpca.mod(mat.sub.t, tol = tolerance, trace = F)


        rm(mat.sub.t)
        gc()
        
        if (verbose){ 
            message("Starting clustering")
        }


        l.mat = rpca.mat$L
        
        rm(rpca.mat)
        gc()
        
        l.mat = t(l.mat)
        rownames(l.mat) = samp.names$nm
        clust.out = hclust(dist(l.mat))
        
        memb = cutree(clust.out, k = number.of.samples)
        memb = setDT(as.data.frame(data.matrix(memb)), keep.rownames = T)
        
        set.seed(12)
        samp.selected = memb[, .SD[sample(1:.N, 1)], by = V1]$rn
        samp.final = normal.table[samp.selected]
        setkey(samp.final, "sample")
    }


    if (use.all){
        if (verbose){ 
            message("Using all samples")
        }


        samp.final = normal.table
        setkey(samp.final, "sample")
    }


    
    message("Balancing pre-decomposition")


    #if (is.human){
        #if (build == "hg19"){
            #par = 2700000
        #} else if (build == "hg38"){
            #par = 2782000
        #} else {
            #stop("provide either hg19 or hg38 build")
        #}
    #}


    if (is.null(PAR.file)){
        message("PAR file not provided, using hg19 default.
If this is not the correct build, please provide a GRange object delineating for corresponding build")
        par.path = system.file("extdata", "PAR_hg19.rds", package = 'dryclean')
        par.gr = readRDS(par.path)
    } else {par.gr = readRDS(PAR.file)}


    message("PAR read")


    samp.final[, file.available := file.exists(normal_cov)]


    message("Checking for existence of files")


    samp.final = samp.final[file.available == TRUE]


    message(paste0(nrow(samp.final), " files present"))
    
    mat.n = pbmclapply(samp.final[, sample], function(nm){
        this.cov = tryCatch(readRDS(samp.final[nm, normal_cov]), error = function(e) NULL)
        chr.prefixed = any('chr' %in% seqnames(this.cov))
        if (!is.null(this.cov)){
            # Handle chromosome prefixed
            all.chr = c(as.character(1:22), "X")
            if (chr.prefixed) {
                all.chr = paste0('chr', c(1:22,'X'))
            }
            #all.chr = c(as.character(1:22), "X")
            ##all.chr = names(which(seqlengths(this.cov) > 5e6))
            this.cov = this.cov %Q% (seqnames %in% all.chr)
            this.cov = this.cov[, field] %>% gr2dt() %>% setnames(., field, "signal.org")
            this.cov = this.cov[!is.na(signal.org),]
            this.cov[, median.idx := .GRP, by = seqnames]
            ##if (is.human){
                ##this.cov[, median.idx := ifelse(seqnames == "X" & start < par, 24, median.idx)]
            ##}
            this.cov$mt = gr.match(dt2gr(this.cov), par.gr)
            this.cov[, median.idx := ifelse(is.na(mt), median.idx, mt+24)]
            ## median.all = this.cov[, .(median.chr = median(signal.org, na.rm = Tp)), by = median.idx]
            ## this.cov = merge(this.cov, median.all, by = "median.idx")
            this.cov[, median.chr := median(signal.org, na.rm = T), by = median.idx]
            this.cov[, signal := ifelse(median.chr == 0, 1, signal.org/median.chr)]
            ## message(nm)
            ## message("this is modified version")
            reads = this.cov[, .(seqnames, signal, median.chr)]
            #reads[, median.chr := median(.SD$signal, na.rm = T), by = seqnames]
            reads[is.na(signal), signal := median.chr]
            min.cov = min(reads[signal > 0]$signal, na.rm = T)
            reads[is.infinite(signal), signal := min.cov]
            reads[signal == 0, signal := min.cov]
            reads[signal < 0, signal := min.cov]
            reads[, signal := log(signal)]
            reads = reads[, .(signal)]
            if (!any(is.infinite(reads$signal))){
                reads = transpose(reads)
                return(reads)
            } 
        } 
    }, mc.cores = num.cores)


    gc()

    mat.bind = rbindlist(mat.n, fill = T)
    mat.bind = na.omit(mat.bind)
    mat.bind.t = transpose(mat.bind)


    rm(mat.bind)
    gc()


    if (verbose){ 
        message("Starting decomposition")
    }


    mat.bind.t = as.matrix(mat.bind.t)
    gc()


    detergent = rrpca.mod(mat.bind.t, trace = F, tol = tolerance)


    rm(mat.bind.t)
    gc()
    
    rsvd.L.burnin = rsvd(detergent$L, k = detergent$k)
    detergent$U.hat = rsvd.L.burnin$u
    detergent$V.hat = t(rsvd.L.burnin$v)
    detergent$sigma.hat = rsvd.L.burnin$d


    if (verbose){ 
        message("Finished making the PON or detergent and saving it to the path provided")
    }


    if (save.pon){
        saveRDS(detergent, paste0(path.to.save, "/detergent.rds"))
    }
    
    return(detergent)

}

rrpca.mod <- function(A, lambda=NULL, maxiter=50, tol=1.0e-3, p=10, q=2, trace=FALSE, rand=TRUE, ...) {
    ##*************************************************************************
    ##***        Author: N. Benjamin Erichson <nbe@st-andrews.ac.uk>        ***
    ##***                              <2016>                               ***
    ##***                       License: BSD 3 clause                       ***
    ##*************************************************************************
    message("This is version 2")
    A <- as.matrix(A)
    m <- nrow(A)
    n <- ncol(A)
    
    rrpcaObj = list(L = NULL,
                    S = NULL,
                    err = NULL)
    
    
    ## Set target rank
    k <- 1
    if(k > min(m,n)) rrpcaObj$k <- min(m,n)
    
    ## Deal with missing values
    is.na(A) <- 0
    
    ## Set lambda, gamma, rho
    if(is.null(lambda)) lambda <- max(m,n)**-0.5
    gamma <- 1.25
    rho <- 1.5
    
    if(rand == TRUE) {
        svdalg = 'rsvd'
    }else {
        svdalg = 'svd'
    }
    
    ## Compute matrix norms
    spectralNorm <- switch(svdalg,
                           svd = norm(A, "2"),
                           rsvd = rsvd(A, k=1, p=10, q=1, nu=0, nv=0)$d,
                           stop("Selected SVD algorithm is not supported!")
                           )
    
    infNorm <- norm( A , "I") / lambda
    dualNorm <- max( spectralNorm , infNorm)
    froNorm <- norm( A , "F")
    
    ## Initalize Lagrange multiplier
    Z <- A / dualNorm
    
    ## Initialize tuning parameter
    mu <- gamma / spectralNorm
    mubar <- mu * 1e7
    mu <- min( mu * rho , mubar )
    muinv <- 1 / mu
    
    ## Init low-rank and sparse matrix
    L = matrix(0, nrow = m, ncol = n)
    S = matrix(0, nrow = m, ncol = n)
    
    niter <- 1
    err <- 1
    while(err > tol && niter <= maxiter) {
        
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ## Update S using soft-threshold
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        epsi = lambda / mu
        temp_S = A - L + Z / mu
        
        S = matrix(0, nrow = m, ncol = n)
        
        idxL <- which(temp_S < -epsi)
        idxH <- which(temp_S > epsi)
        S[idxL] <- temp_S[idxL] + epsi
        S[idxH] <- temp_S[idxH] - epsi
        
        
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ##Singular Value Decomposition
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        R <- A - S + Z / mu
        
        if(svdalg == 'svd') svd_out <- svd(R)
        if(svdalg == 'rsvd') {
            if(k > min(m,n)/5 ) auto_svd = 'svd' else auto_svd = 'rsvd'
            
            svd_out <- switch(auto_svd,
                              svd = svd(R),
                              rsvd = rsvd(R, k=k+10, p=p, q=q))
        }
        
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ## Predict optimal rank and update
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        svp = sum(svd_out$d > 1/mu)
        
        if(svp <= k){
            k = min(svp + 1, n)
        } else {
            k = min(svp + round(0.05 * n), n)
        }
        
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ## Truncate SVD and update L
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ## rrpcaObj$L =  svd_out$u[,1:rrpcaObj$k] %*% diag(svd_out$d[1:rrpcaObj$k] - muinv, nrow=rrpcaObj$k, ncol=rrpcaObj$k)  %*% t(svd_out$v[,1:rrpcaObj$k])
        L =  t(t(svd_out$u[,1:svp, drop=FALSE]) * (svd_out$d[1:svp] - 1/mu)) %*% t(svd_out$v[,1:svp, drop=FALSE])
        
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ## Compute error
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Astar = A - L - S
        Z = Z + Astar * mu
        
        err = norm( Astar , 'F') / froNorm
        rrpcaObj$err <- c(rrpcaObj$err, err)
        
        if(trace==TRUE){
            cat('\n', paste0('Iteration: ', niter ), paste0(' predicted rank = ', svp ), paste0(' target rank k = ', k ),  paste0(' Fro. error = ', rrpcaObj$err[niter] ))
        }
        
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ## Update mu
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        mu = min(mu * rho, mubar);
        muinv = 1 / mu
        
        niter =  niter + 1
        
    }##   End while loop
    rrpcaObj$L <- L
    rrpcaObj$S <- S
    rrpcaObj$k <- k
    
    class(rrpcaObj) <- "rrpca"
    return( rrpcaObj )
    
}


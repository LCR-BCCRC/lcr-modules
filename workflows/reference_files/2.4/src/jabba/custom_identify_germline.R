suppressWarnings(suppressPackageStartupMessages({
    library(GenomicRanges)
    library(data.table)
    library(rsvd)
    library(magrittr)
    library(gUtils)
    library(stats)
    library(MASS)
    library(utils)
    library(pbmcapply)
}))

identify_germline <- function(normal.table.path = NA, signal.thresh = 0.5, pct.thresh = 0.98, verbose = TRUE, save.grm = FALSE, path.to.save = NA, num.cores = 1){


    if (verbose){
        message("Starting the preparation of Panel of Normal samples a.k.a detergent")
    }


    if (is.na(normal.table.path)){
        stop("Need a table with paths to decomposed normal samples to identify germline events")
    }
    
    if (is.na(path.to.save) & save.grm == TRUE){
        stop("Need a path to save identified germline events")
    }


    normal.table = readRDS(normal.table.path)
    setkeyv(normal.table, "sample")
    
    mat.pon = mclapply(normal.table[, sample], function(nm){
        this.cov = tryCatch(readRDS(normal.table[nm, decomposed_cov]), error = function(e) NULL)
        if (!is.null(this.cov)){
            this.cov = gr2dt(this.cov)
            reads = this.cov[, .(foreground.log)] 
            reads = data.table::transpose(reads)
        } else {reads = data.table(NA)}
        return(reads)}, mc.cores = num.cores)
    
    gc()
    
    mat.bind = rbindlist(mat.pon, fill = T)
    ## mat.bind = na.omit(mat.bind)
    mat.bind.t = data.table::transpose(mat.bind)
    
    rm(mat.bind)
    gc()
    
    for(col in names(mat.bind.t)) set(mat.bind.t, i = which(abs(mat.bind.t[[col]]) > signal.thresh), j = col, value = NA) 
    for(col in names(mat.bind.t)) set(mat.bind.t, i = which(!is.na(mat.bind.t[[col]])), j = col, value = 1)
    for(col in names(mat.bind.t)) set(mat.bind.t, i = which(is.na(mat.bind.t[[col]])), j = col, value = 0)


    gc()
    
    mat.bind.t[, black_list_pct := rowSums(.SD)/dim(mat.bind.t)[2]]
    
    mat.bind.t[, germline.status := ifelse(black_list_pct > pct.thresh, FALSE, TRUE)]


    if (nrow(mat.bind.t[germline.status == FALSE]) < 0.5 * nrow(mat.bind.t)){
        warning("More than 50% markers classified as germline, consider adjusting thresholds.")
    }
    
    template = readRDS(normal.table[1, decomposed_cov])
    values(template) <- NULL
    template$germline.status <- mat.bind.t$germline.status


    
    rm(mat.bind.t)
    gc()


    if (verbose){ 
        message("Finished identifying germline markers based pn thresholds provided and saving it to the path provided")
    }


    if (save.grm){
        saveRDS(template, paste0(path.to.save, "/germline.markers.rds"))
    }


    return(template)


}

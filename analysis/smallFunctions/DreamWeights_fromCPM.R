#' @title Estimate Voom Precision Weights Directly From CPM Values, allowing random effects
#'
#' @description Estimate Voom Precision Weights Directly From CPM Values and prepare for Linear Mixed Modelling with \code{dream()}
#'
#' @details
#' This function takes CPM or logCPM values and estimates the precision weights
#' as would be done by providing counts directly to the \code{\link{voomWithDreamWeights}}
#' function.
#' Using this function enables the use of logCPM values which have been
#' normalised using other methods such as Conditional-Quantile Normalisation.
#'
#' The precision weights are returned as part of the output, and these are able
#' to be passed to the function \code{\link{dream}} during linear mixed model fitting.
#' This will ensure that the mean-variance relationship is appropriate for
#' the linear modelling steps as performed by limma.
#'
#'
#' @param cpm Matrix of CPM or logCPM values
#' @param design The model formula containing a random effect
#' @param lib.size Initial library sizes. Must be provided as these are no
#' estimable from CPM values
#' @param isLogCPM logical(1). Indicates whether the data is log2 transformed
#' already. Most commonly (e.g. if using the output of cqn) it will be,
#' @param span Width of the smoothing window used for the lowess mean-variance
#' trend. Expressed as a proportion between 0 and 1.
#' @param BPPARAM parameters for parallel evaluation
#' @param      ... other arguments are passed to \code{lmer}.
#'
#' @return
#' An object of class \code{EList} as would be output by voom.
#' Importantly, there will be no \code{genes} element, although this can be
#' added later.
#' Similarly, the returned \code{targets} element will only contain sample
#' names and library sizes.
#' This can be incorporated with any other metadata as required.
#'
#' Plotting data is always returned, noting the the value \code{sx} has
#' been offset by the library sizes and will be simple logCPM values.
#' As such, the fitted \code{Amean} is also returned in this list element.
#' @import limma
#' @importFrom lme4 VarCorr 
#' @importFrom stats approxfun predict as.formula
#' @export
DreamWeights_fromCPM <- function(cpm, formula, metadata, lib.size=NULL, span=0.5, 
                                 isLogCPM = TRUE, estimateArrayWeight = TRUE, quiet=FALSE, BPPARAM=bpparam(),...){
    
    require(lme4)
    
    formula = as.formula( formula )
    
    out <- list()

    ## Checks taken from voom internals
    n <- nrow(cpm)
    if (n < 2L)
        stop("Need at least two genes to fit a mean-variance trend")
    m <- min(cpm)
    if (is.na(m))
        stop("NA values not allowed")
    if (m < 0 & !isLogCPM)
        stop("Negative CPM values not allowed")
    if (m == 0 & !isLogCPM)
        stop("Please ensure an offset is used for estimation of CPM values.")
    
    ## Make sure we are on the log scale
    if (!isLogCPM)
        cpm <- log2(cpm)
    
    ## Library sizes must be supplied & valid
    if (is.null(lib.size))
        stop("Library sizes must be provded as these cannot estimated from CPM")
    i <- ncol(cpm)
    stopifnot(length(lib.size) == i)
    stopifnot(is.numeric(lib.size) & all(lib.size > 0))
    
    # Fit regression model
    #---------------------
    
    # this function should only be used when random effect exists
     if( missing(metadata) ){
            stop("Must specify argument 'data'\n")
        }
        # fit linear mixed model
        vpList = fitVarPartModel(cpm, formula, metadata, showWarnings=FALSE, ...,fxn = function(fit){
            # extract 
            # 1) sqrt residual variance (i.e. residual standard deviation)
            # 2) fitted values
            list( sd = attr(VarCorr(fit), 'sc'),
                  fitted.values = predict(fit) )
        }, BPPARAM=BPPARAM )
        
        fit = list()
        fit$sigma <- sapply( vpList, function(x) x$sd)	
        fit$df.residual = rep(2, length(fit$sigma)) # check this
        
        # extract fitted values
        fitted.values <- lapply( vpList, function(x) x$fitted.values)
        fitted.values <- do.call("rbind", fitted.values )
        
    if(is.null(fit$Amean)) fit$Amean <- rowMeans(cpm,na.rm=TRUE)
 
    
    # Fit lowess trend to sqrt-standard-deviations by log-count-size
    sx <- fit$Amean+mean(log2(lib.size+1))-log2(1e6)

    # get residual standard deviation
    sy <- sqrt(fit$sigma)
    l <- stats::lowess(sx,sy,f=span)

    suppressWarnings({
        f <- approxfun(l, rule=2)
    })

    fitted.cpm <- 2^fitted.values
    fitted.count <- 1e-6 * t(t(fitted.cpm)*(lib.size+1))
    fitted.logcount <- log2(fitted.count)

    #	Apply trend to individual observations
    w <- 1/f(fitted.logcount)^4
    dim(w) <- dim(fitted.logcount)

    # ## If array weights are provided, the values for `w` would then be used for a
    # ## second round of estimation, and scaled by these weights
    # aw <- c()
    # if (estimateArrayWeight){
    #     ## Use all defaults from this function, providing the precision weights
    #     aw <- arrayWeights(
    #         object = cpm, design = design, weights = w, var.design = NULL,
    #         var.group = NULL, prior.n = 10, method = "auto", maxiter = 50,
    #         tol = 1e-5, trace = FALSE
    #     )
    #     w <- t(aw * t(w))
    # }

    #	Output
    out$E <- cpm
    out$weights <- w
    nm <- colnames(cpm)
    if (is.null(nm)) nm <- seq_len(i)
    out$targets <- data.frame(
        sample = nm,
        lib.size = lib.size
    )
    out$voom.xy <- list(
        x = sx, y = sy, Amean = fit$Amean,
        xlab = "log2( count size + 0.5 )", ylab = "Sqrt ( standard deviation )"
    )
    out$voom.line <- l

    new("EList",out)
  
}

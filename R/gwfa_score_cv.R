gwfa_score_cv <- function(bw, x, dp.locat,k, robust, scores,  elocat=NULL, kernel, adaptive=TRUE, p, theta, longlat, dMat,
                               vars, n.obs = NA,  fm, rotate, oblique.scores=FALSE, timeout, foreach) {

##This function is based on GWmodel::gwpca.cv.
  requireNamespace("GWmodel")
  requireNamespace("psych")
  requireNamespace("foreach")
  #requireNamespace("doMC")
  requireNamespace("doParallel")
  
  data <- x
  if (is(data, "Spatial")) {
    p4s <- proj4string(data)
    dp.locat <- coordinates(data)
  }
  else if (is(data, "data.frame") && (!missing(dMat)))
    data <- data
  else stop("Given data must be a Spatial*DataFrame or data.frame object")
  if (missing(elocat)) {
    ep.given <- FALSE
    elocat <- coordinates(data)
  }
  else {
    ep.given <- T
    if (is(elocat, "Spatial")) {
      espdf <- elocat
      elocat <- coordinates(espdf)
    }
    else if (is.numeric(elocat) && dim(elocat)[2] == 2)
      elocat <- elocat
    else {
      warning("Output loactions are not packed in a Spatial object,and it has to be a two-column numeric vector")
      elocat <- dp.locat
    }
  }
  data <- as(data, "data.frame")
  dp.n <- nrow(data)
  ep.n <- nrow(elocat)
  if (missing(dMat)) {
    DM.given <- F
    DM1.given <- F
    if (dp.n + ep.n <= 10000) {
      dMat <- gw.dist(dp.locat = dp.locat, rp.locat = elocat,
                      p = p, theta = theta, longlat = longlat)
      DM.given <- T
    }
  }
  else {
    DM.given <- T
    DM1.given <- T
    dim.dMat <- dim(dMat)
    if (dim.dMat[1] != dp.n || dim.dMat[2] != ep.n)
      stop("Dimensions of dMat are not correct")
  }
  if (missing(vars))
    stop("Variables input error")
  if (missing(bw) || bw <= 0)
    stop("Bandwidth is not specified incorrectly")
  len.var <- length(vars)
  col.nm <- colnames(data)
  var.idx <- match(vars, col.nm)[!is.na(match(vars, col.nm))]
  if (length(var.idx) == 0)
    stop("Variables input doesn't match with data")
  x <- data[, var.idx]
  x <- as.matrix(x)
  var.nms <- colnames(x)
  var.n <- ncol(x)
  if (len.var > var.n)
    warning("Invalid variables have been specified, please check them again!")


  temp0 <-  tryCatch({ R.utils::withTimeout( 
    fa(x,nfactors = k,n.obs = nrow(x), fm=fm,scores=scores,  rotate=rotate, residuals=T,oblique.scores=oblique.scores),
    timeout=timeout)},
    error=function(e){ NULL})

  if(foreach==TRUE){
    
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    
        cv <- foreach(i= 1:ep.n, .combine = "cbind") %dopar% {
          if (DM.given)
            dist.vi <- dMat[, i]
          else {
            if (ep.given)
              dist.vi <- gw.dist(dp.locat, elocat, focus = i,
                                 p, theta, longlat)
            else dist.vi <- gw.dist(dp.locat, focus = i, p = p,
                                    theta = theta, longlat = longlat)
          }
          wt <- gw.weight(dist.vi, bw, kernel, adaptive)
          use <- wt > 0
          wt[i] <- 0
          wt <- wt[use]
          if (length(wt) <= 5) {
            expr <- paste("Too small bandwidth at location: ",
                          i)
            warning(paste(expr, "and the results can't be given there.",
                          sep = ", "))
            next
          }
          
          temp1 <- wfa(x=data, wt, factors=k, scores=scores, n.obs = length(wt), fm, rotate, oblique.scores=oblique.scores, timeout=timeout)
          
          if(is.null(temp1)){
          out <-  NA
          } else{
          out <-   tryCatch({ sum((temp0$scores[i, ] - temp1$scores[i, ]))**2 },
                     error=function(e){ NA})
          }
          
        }
        
        stopCluster(cl)
        
  } else{
  cv <- c()
  for (i in 1:ep.n) {

    if (DM.given)
      dist.vi <- dMat[, i]
    else {
      if (ep.given)
        dist.vi <- gw.dist(dp.locat, elocat, focus = i,
                           p, theta, longlat)
      else dist.vi <- gw.dist(dp.locat, focus = i, p = p,
                              theta = theta, longlat = longlat)
    }
    wt <- gw.weight(dist.vi, bw, kernel, adaptive)
    use <- wt > 0
    wt[i] <- 0
    wt <- wt[use]
    if (length(wt) <= 5) {
      expr <- paste("Too small bandwidth at location: ",
                    i)
      warning(paste(expr, "and the results can't be given there.",
                    sep = ", "))
      next
    }

    temp1 <- wfa(x=data, wt, factors=k, scores=scores, n.obs = length(wt), fm, rotate, oblique.scores=oblique.scores, timeout=timeout)

    if(is.null(temp1)){
      cv[i] <- NA
    } else{
      cv[i] <- sum((temp0$scores[i, ] - temp1$scores[i, ]))**2
      }
    }
  }
  
  mean(cv %>% as.numeric(), na.rm = TRUE)
  
  }


bw_gwfa <- function(data, vars,k=2, scores, robust=FALSE, kernel, adaptive=TRUE, p=2, theta=0, longlat=FALSE, dMat,
                    n.obs = NA,fm, rotate, type = c("cv_score","cv_uniquenesses", "max_uniquenesses","residual_sum","accumvar_max"),oblique.scores=FALSE, timeout){

  requireNamespace("GWmodel")
  requireNamespace("psych")

  x <- data
  if (is(data, "Spatial")) {
    p4s <- proj4string(data)
    dp.locat <- coordinates(data)
  }
  else if (is(data, "data.frame") && (!missing(dMat)))
    data <- data
  else stop("Given data must be a Spatial*DataFrame or data.frame object")
  data <- as(data, "data.frame")
  dp.n <- nrow(data)
  if (missing(dMat)) {
    DM.given <- F
    if (dp.n <= 5000) {
      dMat <- gw.dist(dp.locat = dp.locat, rp.locat = dp.locat,
                      p = p, theta = theta, longlat = longlat)
      DM.given <- T
    }
  }
  else {
    DM.given <- T
    dim.dMat <- dim(dMat)
    if (dim.dMat[1] != dp.n || dim.dMat[2] != dp.n)
      stop("Dimensions of dMat are not correct")
  }
  if (missing(vars))
    stop("Variables input error")
  col.nm <- colnames(data)
  var.idx <- match(vars, col.nm)[!is.na(match(vars, col.nm))]
  if (length(var.idx) == 0)
    stop("Variables input doesn't match with data")
  data <- data[, var.idx]
  data <- as.matrix(data)
  var.nms <- colnames(data)
  var.n <- ncol(data)
  if (adaptive) {
    upper <- dp.n
    lower <- dp.n/3 ##chenged from 2 to dp.n/3
  }
  else {
    if (DM.given) {
      upper <- range(dMat)[2]
      lower <- upper/5000
    }
    else {
      dMat <- NULL
      if (p == 2) {
        b.box <- bbox(dp.locat)
        upper <- sqrt((b.box[1, 2] - b.box[1, 1])^2 +
                        (b.box[2, 2] - b.box[2, 1])^2)
        lower <- upper/5000
      }
      else {
        upper <- 0
        for (i in 1:dp.n) {
          dist.vi <- gw.dist(dp.locat = dp.locat, focus = i,
                             p = p, theta = theta, longlat = longlat)
          upper <- max(upper, range(dist.vi)[2])
        }
        lower <- upper/5000
      }
    }
  }
  bw <- NA
  if(type=="cv_uniquenesses"){
    bw <- gold(gwfa.cv_uniquenesses.calc, lower, upper, adapt.bw = adaptive, x,
               dp.locat, k, elocat=NULL, robust, kernel, adaptive, p, theta, longlat,
               dMat, vars,fm=fm,rotate=rotate,scores=scores,oblique.scores=oblique.scores, timeout=timeout)
  } else if (type=="max_uniquenesses"){
    bw <- gold(gwfa_uniquenesses_sum, lower, upper, adapt.bw = adaptive, x,
               dp.locat, k, robust, kernel, adaptive, p, theta, longlat,
               dMat, vars,fm=fm,rotate=rotate,scores=scores,oblique.scores=oblique.scores, timeout=timeout)
  } else if (type=="cv_score"){
    bw <- gold(gwfa_score_cv, lower, upper, adapt.bw = adaptive, x,
               dp.locat, k, elocat = NULL, robust, kernel, adaptive, p, theta, longlat,
               dMat, vars,fm=fm,rotate=rotate,scores=scores,oblique.scores=oblique.scores,timeout=timeout)
  } else if (type=="residual_sum"){
    bw <- gold(gwfa_residual_sum, lower, upper, adapt.bw = adaptive, x,
               dp.locat, k,  robust, kernel, adaptive, p, theta, longlat,
               dMat, vars,fm=fm,rotate=rotate,scores=scores,oblique.scores=oblique.scores,timeout=timeout)
  } else if (type=="accumvar_max"){
    bw <- gold(gwfa.Accumvar_max.calc, lower, upper, adapt.bw = adaptive, x,
               dp.locat, k,  robust, kernel, adaptive, p, theta, longlat,
               dMat, vars,fm=fm,rotate=rotate,scores=scores,oblique.scores=oblique.scores,timeout=timeout)
  } else {bw <- NA }

  bw
}

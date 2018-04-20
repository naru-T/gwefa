gwefa <- function(data,elocat, vars,bw,k=2, kernel, adaptive=TRUE, p=2, theta=0, longlat=FALSE, dMat=NULL,
                  n.obs = NA, n.iter=1, rotate="oblimin", scores="regression",
                  residuals=FALSE, SMC=TRUE, covar=FALSE,missing=FALSE,impute="median",
                  min.err = 0.001,  max.iter = 50,symmetric=TRUE, warnings=TRUE, fm="minres",
                  alpha=.1,pr=.05,oblique.scores=FALSE,np.obs=NULL,use="pairwise",cor="cor",
                  correct=.5,weight=NULL,...) {


##This function is based on GWmodel::gwpca and psych::fa
  requireNamespace("GWmodel")
  requireNamespace("psych")

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

  load <- array(NA,c(ep.n,var.n,k))
  s <- matrix(NA,ep.n,k)
  u <- matrix(NA,ep.n,var.n)
  ld <- matrix(NA,ep.n,k)
  ss <- matrix(NA,ep.n,k)
  cor.mt <- matrix(NA,ep.n,sum(seq(1,(k-1),1)))
  rmsea <- rep(NA, nrow(elocat))
  score.all <-  matrix(NA,ep.n,k)

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
    wt <- wt[use]
    if (length(wt) <= 5) {
      expr <- paste("Too small bandwidth at location: ",
                    i)
      warning(paste(expr, "and the results can't be given there.",
                    sep = ", "))
      next
    }

    temp <- wefa(x=data[use, ], wt, factors=k, scores=scores, n.obs, fm=fm, rotate=rotate)

    load[i,,] <- matrix(temp$loadings, ncol=k, nrow=var.n)
    #score
    score.all[use, ] <- temp$scores
    s[i,] <- score.all[i,]
    u[i,] <- temp$uniquenesses
    ld[i,]<- temp$Vaccounted[1,]
    ss[i,] <- temp$Vaccounted[2,]
    cor.mt[i,] <- temp$r.scores[upper.tri(temp$r.scores)]
    rmsea[i] <- temp$RMSEA[1]
  }
  dimnames(load)[[3]] <- paste0("Factor",1:k)
  list(
    loadings=load,
    score=s,
    uniquenesses=u,
    loading_score=ld,
    ss=ss,
    cor.mt=cor.mt,
    rmsea=rmsea)
}

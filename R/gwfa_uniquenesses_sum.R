gwfa_uniquenesses_sum <- function(bw, x ,dp.locat, k, robust, kernel, adaptive, p, theta, longlat, dMat,vars, n.obs, fm, rotate,scores, oblique.scores, timeout, foreach){

  ans <- gwfa(data=x,elocat=dp.locat, vars,bw,k, kernel, adaptive, p=2, theta=0, longlat, dMat, n.obs, fm, rotate,scores,oblique.scores=oblique.scores, timeout=timeout, foreach=foreach)

  return(-sum(ans$uniquenesses^2))

}

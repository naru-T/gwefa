gwfa_residual_sum <- function(bw, x ,dp.locat, k, robust, kernel, adaptive, p, theta, longlat, dMat,vars, n.obs, fm, rotate,scores,oblique.scores=FALSE, timeout){

  ans <- gwfa(data=x,elocat=dp.locat, vars,bw,k, kernel, adaptive, p=2, theta=0, longlat, dMat, n.obs, fm, rotate,scores, oblique.scores=oblique.scores, timeout=timeout)

  return(sum(ans$residuals_sqsum) + sum(ans$uniquenesses^2) )

}

gwfa.Accumvar_max.calc <- function(bw, x ,dp.locat, k, robust, kernel, adaptive, p, theta, longlat, dMat,vars, n.obs, fm, rotate,scores, oblique.scores=FALSE, timeout){

  requireNamespace("dplyr")

  temp0 <- tryCatch({ R.utils::withTimeout( 
    fa(data.frame(x) %>% data.frame() %>%
         dplyr::select(vars) ,nfactors=k, fm=fm,rotate=rotate, scores=scores, residuals=T, oblique.scores=oblique.scores),
    timeout=timeout)},
    error=function(e){ NULL})
  temp <- gwfa(data=x,elocat=dp.locat, vars,bw,k, kernel, adaptive, p=2, theta=0, longlat, dMat, n.obs, fm, rotate,scores, oblique.scores=oblique.scores, timeout=timeout)

  accumv <- (1 - min(rowSums(temp$ss)))
  #accumv <- (1- mean(rowSums(temp$ss)))
  #accumv <- (1- median(rowSums(temp$ss)))
  #accumv <- (temp0$Vaccounted[2,] - rowSums(temp$ss))^2
  if(is.na(accumv)){
    return(NA)
  } else{return(accumv)}
}

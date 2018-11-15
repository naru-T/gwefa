gwfa.score_dif.calc <- function(bw, x ,dp.locat, k, robust, kernel, adaptive, p, theta, longlat, dMat,vars, n.obs, fm, rotate,scores,timeout){
  
  temp0 <- fa(x %>% data.frame() %>% 
                dplyr::select(vars) ,nfactors=k, fm=fm,rotate=rotate, scores=scores, residuals=T,timeout=timeout)
  temp1 <- gwfa(data=x,elocat=dp.locat, vars,bw,k, kernel, adaptive, p=2, theta=0, longlat, dMat, n.obs, fm, rotate,scores,timeout=timeout)
  
  sum((temp0$scores - temp1$score)**2)
    
}

gwfa.score_dif.calc <- function(bw, x ,dp.locat, k, robust, kernel, adaptive, p, theta, longlat, dMat,vars, n.obs, fm, rotate,scores){
  
  temp0 <- fa(x %>% data.frame() %>% 
                dplyr::select(vars) ,nfactors=k, fm=fm,rotate=rotate, scores=scores, residuals=T)
  temp1 <- gwfa(data=x,elocat=dp.locat, vars,bw,k, kernel, adaptive, p=2, theta=0, longlat, dMat, n.obs, fm, rotate,scores)
  
  sum((temp0$scores - temp1$score)**2)
    
}

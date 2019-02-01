wfa <- function(x, wt, factors, n.obs = NA, scores="regression", fm = "minres", rotate = "Promax", control = NULL, oblique.scores, timeout,...) {
  ##This function is based on psych::fa
  requireNamespace("psych")

  local.center <- function(x,wt){
    sweep(x,2,colSums(sweep(x,1,wt,'*'))/sum(wt))}
  localx <- sweep(local.center(x,wt),1,sqrt(wt),'*')
  
  out <- tryCatch({ R.utils::withTimeout( 
    fa(localx,nfactors = factors,n.obs = n.obs, fm = fm, scores=scores, rotate = rotate, residuals=T, oblique.scores=oblique.scores),
    timeout=timeout)},
    error=function(e){ NULL})
  return(out)
}

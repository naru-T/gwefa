wefa <- function(x, wt, factors, n.obs = NA, scores="regression", fm = "minres", rotate = "Promax", control = NULL, ...) {
  ##This function is based on psych::fa
  requireNamespace("psych")

  local.center <- function(x,wt){
    sweep(x,2,colSums(sweep(x,1,wt,'*'))/sum(wt))}
  localx <- sweep(local.center(x,wt),1,sqrt(wt),'*')
  fa(localx,nfactors = factors,fm = fm, scores=scores, rotate = rotate, residuals=T)

}

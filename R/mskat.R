#' Fit a null multi-trait regression model
#'
#' Compute a null multivariate regression model
#' @param Ys  matrix of continuous outcomes with variables in column
#' @param X   matrix of covariates to be adjusted
#' @export
MSKAT.cnull <- function(Ys,X=NULL){
  N = dim(Ys)[1]
  ## null model fitting
  if(is.null(X)){
    px = 1
    Ux = matrix(1/sqrt(N), N,1)
  } else{
    px = dim(X)[2]+1
    Ux = svd(cbind(1,X), nv=0)$u
  }
  U0 = Ys - Ux%*%(t(Ux)%*%Ys)
  Sig = cor(U0); iSig = solve(Sig)
  a1 = svd(Sig,nv=0); ihSig = a1$u%*%diag(1/sqrt(a1$d))%*%t(a1$u)
  V0 = colSums(U0^2)/(N-px)
  return(list(Ux=Ux,U0=U0,V0=V0,Sig=Sig,iSig=iSig,ihSig=ihSig,svd=a1, Ys=Ys,px=px) )
}

#' Multi-trait sequence kernel association test of variant set
#'
#' Compute the variant set association test statistic and p-value
#' @param obj  fitted null model from MSKAT.cnull
#' @param G    genotype matrix with SNPs in column
#' @param W    variant weight; use Beta dist weight if NULL
#' @param W.beta   Beta distribution parameters for the variant weight
#' @return
#' \describe{
#'   \item{p.value}{ p-values for the two kernel statistics (Q and Q'; see ref)}
#'   \item{Q}{ test statistics (Q and Q'; see ref)}
#' }
#' @export
#' @references
#' Wu,B. and Pankow,J.S. (2016) Sequence kernel association test of multiple continuous phenotypes. Genetic Epidemiology, 40(2), 91-100.
#' @examples
#' library(CompQuadForm)
#' Y = matrix(rnorm(2000), 1000,2)
#' X = matrix(rnorm(2000), 1000,2)
#' G = matrix(rbinom(10000,2,0.02), 1000,10)
#' MSKAT(MSKAT.cnull(Y,X), G, W.beta=c(1,25))
#' ## time-consuming; not run
#' ## X.list = vector('list',2); for(i in 1:2) X.list[[i]] = rbind(1, t(X))
#' ## MKMR(t(Y),X.list, Gm=G, W.beta=c(1,25))
MSKAT <- function(obj, G, W=NULL, W.beta=c(1,25)){
  n = dim(G)[1];  m = dim(G)[2]; K = dim(obj$Ys)[2]
  if(is.null(W)){
    W = dbeta(colMeans(G)/2, W.beta[1],W.beta[2])
    W = W/sum(W)*m
  } else{
    if(length(W)<=1) W = rep(1,m)
  }
  ## KAT statistic
  Ux = obj$Ux; U0 = obj$U0; V0 = obj$V0
  S = t(t(U0)%*%G/sqrt(V0)) ## m,K
  Ge = G - Ux%*%(t(Ux)%*%G); Ge2 = t(Ge)%*%Ge
  R = t(Ge2*W)*W
  Lam = eigen(R, sym=TRUE, only.val=TRUE)$val
  ##
  ihSig = obj$ihSig
  Z = W*S%*%ihSig
  Q = sum(Z^2)
  pval = KATpval(Q,rep(Lam,K))
  ##
  Z = W*S
  Lam1 = eigen(obj$Sig,sym=TRUE,only.val=TRUE)$val 
  Qp = sum(Z^2)
  pvalp = KATpval(Qp, c(outer(Lam,Lam1)))
  return( list(p.value=c(pval, pvalp), Q=c(Q,Qp)) )
}


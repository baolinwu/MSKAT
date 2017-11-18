### saddlepoint approx: modified from Lumley survey package.
saddle = function(x,lambda){
  d = max(lambda)
  lambda = lambda/d
  x = x/d
  k0 = function(zeta) -sum(log(1-2*zeta*lambda))/2
  kprime0 = function(zeta) sapply(zeta, function(zz) sum(lambda/(1-2*zz*lambda)))
  kpprime0 = function(zeta) 2*sum(lambda^2/(1-2*zeta*lambda)^2)
  n = length(lambda)
  if (any(lambda < 0)) {
    lmin = max(1/(2 * lambda[lambda < 0])) * 0.99999
  } else if (x>sum(lambda)){
    lmin = -0.01
  } else {
    lmin = -length(lambda)/(2*x)
  }
  lmax = min(1/(2*lambda[lambda>0]))*0.99999
  hatzeta = uniroot(function(zeta) kprime0(zeta) - x, lower = lmin, upper = lmax, tol = 1e-08)$root
  w = sign(hatzeta)*sqrt(2*(hatzeta*x-k0(hatzeta)))
  v = hatzeta*sqrt(kpprime0(hatzeta))
  if(abs(hatzeta)<1e-4){
    return(NA)
  } else{
    return( pnorm(w+log(v/w)/w, lower.tail=FALSE) )
  }
}
Sadd.pval = function(Q.all,lambda){
  sad = rep(1,length(Q.all))
  if(sum(Q.all>0)>0){
    sad[Q.all>0] = sapply(Q.all[Q.all>0],saddle,lambda=lambda)
  }
  id = which(is.na(sad))
  if(length(id)>0){
    sad[id] = Liu.pval(Q.all[id], lambda)
  }
  return(sad)
}
### modified Liu method from Lee SKAT-O paper
Liu.pval = function(Q, lambda){
  c1 = rep(0,4); for(i in 1:4){ c1[i] = sum(lambda^i) }
  muQ = c1[1];  sigmaQ = sqrt(2 *c1[2])
  s1 = c1[3]/c1[2]^(3/2);  s2 = c1[4]/c1[2]^2
  if(s1^2 > s2){
    a = 1/(s1 - sqrt(s1^2 - s2));  d = s1 *a^3 - a^2;  l = a^2 - 2*d
  } else {
    l = 1/s2;  a = sqrt(l);  d = 0
  }
  muX = l+d;  sigmaX = sqrt(2)*a
  
  Q.Norm = (Q - muQ)/sigmaQ
  Q.Norm1 = Q.Norm*sigmaX + muX
  pchisq(Q.Norm1, df = l,ncp=d, lower.tail=FALSE)
}
#' Compute the tail probability of 1-DF chi-square mixtures
#'
#' Use Davies' method; if fail, switch to the Saddlepoint approx 
#' @param Q.all  test statistics
#' @param lambda  mixing coefficients
#' @param acc  error bound
#' @param lim  maximum number of integration terms
#' @export
KAT.pval <- function(Q.all, lambda, acc=1e-9,lim=1e6){
  pval = rep(0, length(Q.all))
  i1 = which(is.finite(Q.all))
  for(i in i1){
    tmp = davies(Q.all[i],lambda,acc=acc,lim=lim); pval[i] = tmp$Qq
    if((tmp$ifault>0)|(pval[i]<=0)|(pval[i]>=1)) pval[i] = Sadd.pval(Q.all[i],lambda)
  }
  return(pval)
}

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
  Lam = svd(R, nu=0,nv=0)$d
  Lam = Lam[Lam>0]
  ##
  ihSig = obj$ihSig
  Z = W*S%*%ihSig
  ## Zv = kronecker(R,diag(K))
  Q = sum(Z^2)
  pval = KAT.pval(Q,rep(Lam,K))
  ##
  Z = W*S
  ## Zv = kronecker(R, obj$Sig)
  ## Lam1 = svd(obj$Sig,nu=0,nv=0)$d; 
  Lam1 = eigen(obj$Sig,sym=TRUE,only.val=TRUE)$val 
  Lam1 = Lam1[Lam1>0]
  Qp = sum(Z^2)
  pvalp = KAT.pval(Qp, c(outer(Lam,Lam1)))
  return( list(p.value=c(pval, pvalp), Q=c(Q,Qp)) )
}


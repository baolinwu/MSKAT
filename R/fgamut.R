#' Fast implementation of GAMuT (Gene Association with Multiple Traits)
#'
#' Return the multi-trait variant set association test p-values from the Broadaway et. al (2016) approach.
#'  Both tests from linear kernel and projection matrix are computed.
#' This is a much faster implementation following the approach of Wu and Pankow (2016) using the multiple
#' linear regression based score test. We also provide modified (and more accurate) p-values following the approach of
#' Wu et al. (2016). 
#' @param Y  outcome matrix. sample in rows.
#' @param G   genotype matrix. sample in rows.
#' @param X  covariate matrix to be adjusted. sample in rows.
#' @param W    variant weight; use Beta dist weight if NULL
#' @param W.beta   Beta distribution parameters for the variant weight. Default to Beta(1,25).
#' @return a p-value vector p.value with four components
#' \describe{
#'   \item{proj}{ p-value for the GAMuT with projection matrix kernel }
#'   \item{Cproj}{ corrected p-value for the GAMuT with projection matrix kernel }
#'   \item{linear}{ p-value for the GAMuT with linear kernel }
#'   \item{Clinear}{ corrected p-value for the GAMuT with linear kernel }
#' }
#' @export
#' @references
#' Broadaway, K.A., Cutler, D.J., Duncan, R., Moore, J.L., Ware, E.B., Jhun, M.A., Bielak, L.F., Zhao, W., Smith, J.A., Peyser, P.A., Kardia, S.L.R., Ghosh, D., Epstein, M.P. (2016). A Statistical Approach for Testing Cross-Phenotype Effects of Rare Variants. \emph{The American Journal of Human Genetics}, 98, 525--540.
#'
#' Wu,B. and Pankow,J.S. (2016) Sequence kernel association test of multiple continuous phenotypes. \emph{Genetic Epidemiology}, 40(2), 91-100.
#' 
#' Wu,B., Guan,W., Pankow,J.S. (2016) On efficient and accurate calculation of significance p-values for sequence kernel association test of variant set. \emph{Annals of human genetics}, 80(2), 123-135.
#' 
#' Wu,B. and Qin,H. (2018). On testing cross-phenotype effects of rare variants. \emph{AHG}, under revision. \url{http://www.umn.edu/~baolin/research/gamut}
fgamut <- function(Y,G, X=NULL, W=NULL, W.beta=c(1,25)){
  if(class(Y)!='matrix') Y = as.matrix(Y)
  if(class(G)!='matrix') G = as.matrix(G)
  n = dim(Y)[1]; K = dim(Y)[2]; m = dim(G)[2]
  if(is.null(W)){
    W = dbeta(colMeans(G)/2, W.beta[1],W.beta[2])
    W = W/sum(W)*m
  } else{
    if(length(W)<=1) W = rep(1,m)
  }
  G = scale(t(t(G)*W), TRUE, FALSE)
  lamg = zapsmall( svd(G,nu=0,nv=0)$d^2 )
  lamg = lamg[lamg>0]
  ## residuals calc
  if(is.null(X)){
    resid = Y
    p = 1
  } else{
    if(class(X)!='matrix') X = as.matrix(X)
    p = dim(X)[2]+1
    Ux = svd(cbind(1,X),nv=0)$u
    resid = Y - Ux%*%(t(Ux)%*%Y)
  }
  P0 = scale(resid, TRUE,TRUE)
  a0 = svd(P0,nv=0)
  Uy = a0$u
  lamy = a0$d^2/n
  ## linear kernel
  S = t(G)%*%P0 ## m,K
  Ql = sum(S^2)
  Laml = c(outer(lamy,lamg))
  pvall = davies(Ql,Laml)$Qq
  pval.l = KAT.pval(Ql,Laml)
  ## proj matrix
  Sp = t(G)%*%Uy
  Qp = sum(Sp^2)
  Lamp = rep(lamg,K)/n
  pval.p = KAT.pval(Qp,Lamp)
  pvalp = davies(Qp,Lamp)$Qq
  p.value = c(pvalp,pval.p, pvall,pval.l)
  names(p.value) = c('proj','Cproj','linear','Clinear')
  return(p.value)
}

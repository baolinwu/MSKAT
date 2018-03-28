#' GWAS summary data based SNP-set association test across multiple phenotypes
#'
#' Efficient multi-trait SNP-set association tests: sum of squared test (S2T), squared burden test (SBT), their adaptive test:
#' Fisher combination test (FCT), and Tippett combination test (TCT).
#' 
#' @param  Z  M by K matrix of summary Z-statistics for M variants across K traits
#' @param  Sig estimated trait correlation matrix (K by K)
#' @param  R  variant LD correlation matrix (M by M)
#' @return p-value vector: FCT, SBT, S2T and TCT
#' @export
#' @references
#' Guo,B. and Wu,B. (2018) GWAS summary data based SNP-set association test across multiple phenotypes to identify novel genetic variants. tech rep.
#' @examples
#' K = 4; M = 20
#' R = cor(matrix(rnorm(500*M),500,M)*sqrt(0.8)+rnorm(500)*sqrt(0.2))
#' Y = matrix(rnorm(100*K), 100,K)*sqrt(0.8)+rnorm(100)*sqrt(0.2)
#' Y[,1] = -Y[,1]; Sig = cor(Y)
#' Z = matrix(rnorm(K*M),M,K)
#' js = sample(1:M, size=round(M*0.4))
#' Z[js,] = Z[js,] + 1.5
#' msatz(Z,Sig,R)
#' Z = matrix(rnorm(K*M),M,K)
#' ij = cbind(sample(1:M, size=20, rep=TRUE), sample(1:K, size=20, rep=TRUE))
#' Z[ij] = Z[ij] + rnorm(20)*1.5
#' msatz(Z,Sig,R)
msatz <- function(Z, Sig, R){
  if(class(Z)!='matrix') Z = as.matrix(Z)
  M = dim(Z)[1]; K = dim(Z)[2]
  lamR = eigen(R,sym=TRUE,only.val=TRUE)$val
  lamS = eigen(Sig,sym=TRUE,only.val=TRUE)$val
  ## S2T
  pvalv = KATpval(sum(Z^2), outer(lamR,lamS))
  R1 = rowSums(R); S1 = rowSums(Sig)
  R2 = sum(R1); S2 = sum(S1)
  ## SBT
  B = colSums(Z); Qb = sum(B^2)
  pvalb = KATpval(Qb/R2, lamS)
  ## E
  E = Z - outer(R1,B)/R2
  Vr = R - outer(R1,R1)/R2
  lamr = eigen(Vr, sym=TRUE,only.val=TRUE)$val
  pvale = KATpval(sum(E^2), outer(lamS,lamr) )
  ## Fisher comb
  Qb = -2*log(pvalb); Qe = -2*log(pvale)
  fval = pchisq(Qb+Qe, 4, lower=FALSE)
  ## Tippet (minP)
  tval = 1-(1-min(pvale,pvalb))^2
  return( list(p.value=c(FCT=fval, SBT=pvalb, S2T=pvalv, TCT=tval)) )
}

#' GWAS summary data based phenome-wide SNP-set association test
#'
#' Efficient and powerful summary data based phenome-wide SNP-set association tests: quadratic test (S2T), burden test (SBT), their adaptive test:
#' Fisher combination test (FCT), and Tippett combination test (TCT).
#' 
#' @param  Z  M by K matrix of summary Z-statistics for M variants across K traits
#' @param  Sig estimated trait correlation matrix (K by K)
#' @param  R  variant LD correlation matrix (M by M)
#' @return p-value vector: FCT, SBT, S2T, TCT, and S2E (residual effects excluding SBT)
#' @export
#' @references
#' Guo,B., Bishop,J., Liu,N. and Wu,B. (2018) Phenome-wide SNP-set association test based on GWAS summary data to identify novel disease-gene association. tech rep.
#' @examples
#' K = 4; M = 20
#' R = cor(matrix(rnorm(500*M),500,M)*sqrt(0.8)+rnorm(500)*sqrt(0.2))
#' Y = matrix(rnorm(100*K), 100,K)*sqrt(0.8)+rnorm(100)*sqrt(0.2); Y[,1] = -Y[,1]
#' Sig = cor(Y)
#' Rc = t(chol(R)); Sc = chol(Sig)
#' Z0 = Rc%*%matrix(rnorm(K*M),M,K)%*%Sc
#' PheSATS(Z0,Sig,R)
#' js = sample(1:M, size=round(M*0.4))
#' Z1 = Z0; Z1[js,] = Z1[js,] + 1.5
#' PheSATS(Z1,Sig,R)
#' ij = cbind(sample(1:M, size=20, rep=TRUE), sample(1:K, size=20, rep=TRUE))
#' Z2 = Z0; Z2[ij] = Z2[ij] + rnorm(20)*1.5
#' PheSATS(Z2,Sig,R)
PheSATS <- function(Z, Sig, R){
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
  ## Tippett
  tval = 1-(1-min(pvale,pvalb))^2
  return( list(p.value=c(FCT=fval, SBT=pvalb, S2T=pvalv, TCT=tval, S2E=pvale)) )
}

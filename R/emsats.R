#' PCA based adapitve SNP-set association test of multiple traits using GWAS summary statistics
#'
#' We use the LD score regression (Bulik-Sullivan et al.) to accurately estimate the trait correlation,
#' which is then used to construct the multi-trait association tests of multiple variants in a gene or pathway: PC based test (ET), 
#' variance components test (VT), and their adaptive test (AT)
#' 
#' @param  Z  M by K matrix of summary Z-statistics for M variants across K traits
#' @param  Sig estimated trait correlation matrix (K by K)
#' @param  R  variant LD correlation matrix (M by M)
#' @param  rho  weight assigned to the ET
#' @return
#' \describe{
#'   \item{p.value}{ test p-values for: AT, VC, ET }
#'   \item{pval}{ vector of all p-values }
#'   \item{rho.est}{ the optimal rho weight }
#' }
#' @export
#' @references
#' Guo,B. and Wu,B. (2017) PCA based adaptive multi-trait SNP-set association test using the GWAS summary statistics. tech rep.
#' @examples
#' K = 4; M = 20
#' R = cor(matrix(rnorm(500*M),500,M)*sqrt(0.8)+rnorm(500)*sqrt(0.2))
#' Y = matrix(rnorm(100*K), 100,K)*sqrt(0.8)+rnorm(100)*sqrt(0.2)
#' Y[,1] = -Y[,1]; Sig = cor(Y)
#' Z = matrix(rnorm(K*M),M,K)
#' js = sample(1:M, size=round(M*0.4))
#' Z[js,] = Z[js,] + 1.5
#' EMSATS(Z,Sig,R)
#' Z = matrix(rnorm(K*M),M,K)
#' ij = cbind(sample(1:M, size=20, rep=TRUE), sample(1:K, size=20, rep=TRUE))
#' Z[ij] = Z[ij] + rnorm(20)*1.5
#' EMSATS(Z,Sig,R)
EMSATS <- function(Z,Sig, R, rho=c(-2,-1,0,0.5,1)){
  if(class(Z)!='matrix') Z = as.matrix(Z)
  M = dim(Z)[1]; K = dim(Z)[2]
  a0 = eigen(Sig, sym=TRUE)
  u1d = a0$vec[,1]/sqrt(a0$val[1])
  iSig = solve(Sig)
  lamR = eigen(R, sym=TRUE, only.val=TRUE)$val
  ##
  Q1 = sum(Z%*%iSig*Z)
  Qb = sum( colSums(u1d*t(Z))^2 )
  Qs = (1-rho)*Q1 + rho*Qb
  ## 
  L = length(rho); pval = rep(1,L)
  LAM = vector('list', L)
  lam1 = rep(lamR, K-1)
  for(j in 1:L){
    LAM[[j]] = c( (1-rho[j])*lam1, lamR)
    pval[j] = KATpval(Qs[j], LAM[[j]])
  }
  minP = min(pval)
  ##
  qval = rep(1,L)
  for(j in 1:L){
    qval[j] = KATqval(minP, LAM[[j]])
  }
  B=1e3; q1 = qval[L]
  q1x = seq(0, q1, length=B)
  dx = sapply(q1x, function(x) min((qval[-L]-x)/(1-rho[-L])))
  p1 = KATpval(dx, lam1);  p2 = KATpval(q1x, lamR)
  p.val = minP - sum((p1[-1]+p1[-B])/2*diff(p2))
  return( list(p.value=c(AT=pvalo,VC=pval[rho==0],ET=pval[rho==1]), pval=pval, rho.est=rho[which.min(pval)]) )
}


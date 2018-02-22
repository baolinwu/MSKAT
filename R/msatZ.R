#' Multi-trait SNP-set association test using GWAS summary data
#'
#' Efficient multi-trait SNP-set association tests: squared sum test (SST), sum test (ST), their adaptive test (AT)
#' 
#' @param  Z  M by K matrix of summary Z-statistics for M variants across K traits
#' @param  Sig estimated trait correlation matrix (K by K)
#' @param  R  variant LD correlation matrix (M by M)
#' @param  rho  weight assigned to the ST
#' @return
#' \describe{
#'   \item{p.value}{ test p-values for: AT, SST, ST }
#'   \item{pval}{ vector of all p-values }
#'   \item{rho.est}{ the optimal rho weight }
#' }
#' @export
#' @references
#' Guo,B. and Wu,B. (2017) GWAS summary data based multi-trait SNP-set association test to identify novel genetic variants. tech rep.
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
msatz <- function(Z,Sig, R, rho=0:10/10){
  if(class(Z)!='matrix') Z = as.matrix(Z)
  M = dim(Z)[1]; K = dim(Z)[2]
  ##
  Zw = as.vector(Z); Rw = kronecker(Sig,R)
  eR = eigen(Rw,sym=TRUE);  lamR = abs(eR$val); eta = colSums(eR$vec)*sqrt(lamR)
  R1 = sum(eta^2); R2 = sum(eta^2*lamR)
  c2 = outer(eta,eta)
  Lamq = eigen(diag(lamR) - R2/R1^2*c2, symmetric=TRUE,only.values=TRUE)$val
  ##
  ## lamR = eigen(R, sym=TRUE, only.val=TRUE)$val
  ## lamS = eigen(Sig, sym=TRUE, only.val=TRUE)$val
  ## ST
  Qs = sum(Zw)^2
  pvals = pchisq(Qs/sum(Sig)/sum(R), 1,lower=FALSE)
  ## SST
  Q2 = sum(Zw^2); 
  pval2 = KATpval(Q2,lamR)
  ## AT
  L = length(rho)
  if(L<=2){
    return(list(p.value=c(A=NULL, S2=pval2, S=pvals), pval=c(pval2,pvals)) )
  } 
  L1 = L-1; rho1 = rho[-L]
  Qw = (1-rho)*Q2 + rho*Qs
  pval = rep(1, L)
  pval[1] = pval2; pval[L] = pvals
  Lamk = vector('list', L)
  Lamk[[L]] = R1;  Lamk[[1]] = lamR
  for(k in 2:L1){
    mk = rho[k]*c2;  diag(mk) = diag(mk) + (1-rho[k])*lamR
    aak = zapsmall( abs(eigen(mk, sym=TRUE, only.val=TRUE)$val) )
    Lamk[[k]] = aak[aak>0]
    pval[k] = KATpval(Qw[k],Lamk[[k]])
  }
  minP = min(pval)
  L = length(rho)
  qval = rep(0,L1)
  for (k in 1:L1) qval[k] = Liu0.qval(minP, Lamk[[k]])
  q1 = qchisq(minP,1,lower=FALSE)
  tauk = (1-rho1)*R2/R1 + rho1*R1
  katint = function(xpar){
    eta1 = sapply(xpar, function(eta0) min((qval-tauk*eta0)/(1-rho1)))
    KATpval(eta1,Lamq)*dchisq(xpar,1)
  }
  prec = 1e-4
  p.value = try({ minP + integrate(katint, 0,q1,  subdivisions=1e3,abs.tol=minP*prec)$val }, silent=TRUE)
  while(class(p.value)=='try-error'){
    prec = prec*2
    p.value = try({ minP + integrate(katint, 0,q1, abs.tol=minP*prec)$val }, silent=TRUE)
  }
  p.value = min(p.value,minP*L)
  return(list(p.value=c(A=p.value, S2=pval2, S=pvals), pval=pval, rho.est=rho[which.min(pval)]) )
}

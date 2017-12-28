#' Multi-trait SNP-set association test using GWAS summary data
#'
#' Efficient multi-trait SNP-set association tests: variance components test (VC), burden type test (BT), adaptive test (AT)
#' 
#' @param  Z  M by K matrix of summary Z-statistics for M variants across K traits
#' @param  Sig estimated trait correlation matrix (K by K)
#' @param  R  variant LD correlation matrix (M by M)
#' @param  rho  weight assigned to the BT
#' @return
#' \describe{
#'   \item{p.value}{ test p-values for: AT, VC, BT }
#'   \item{pval}{ vector of all p-values }
#'   \item{rho.est}{ the optimal rho weight }
#' }
#' @export
#' @references
#' Guo,B. and Wu,B. (2017) Powerful and efficient SNP-set association tests across multiple phenotypes using GWAS summary data. tech rep.
#' @examples
#' K = 4; M = 20
#' R = cor(matrix(rnorm(500*M),500,M)*sqrt(0.8)+rnorm(500)*sqrt(0.2))
#' Y = matrix(rnorm(100*K), 100,K)*sqrt(0.8)+rnorm(100)*sqrt(0.2)
#' Y[,1] = -Y[,1]; Sig = cor(Y)
#' Z = matrix(rnorm(K*M),M,K)
#' js = sample(1:M, size=round(M*0.4))
#' Z[js,] = Z[js,] + 1.5
#' MSATS(Z,Sig,R)
#' Z = matrix(rnorm(K*M),M,K)
#' ij = cbind(sample(1:M, size=20, rep=TRUE), sample(1:K, size=20, rep=TRUE))
#' Z[ij] = Z[ij] + rnorm(20)*1.5
#' MSATS(Z,Sig,R)
MSATS <- function(Z,Sig, R, rho=c((0:5/10)^2,0.5,1)){
  if(class(Z)!='matrix') Z = as.matrix(Z)
  M = dim(Z)[1]; K = dim(Z)[2]
  ##
  Rs = rowSums(R); R1 = sum(Rs); R2 = sum(Rs^2); R3 = sum(Rs*colSums(R*Rs))
  RJ2 = outer(Rs,Rs,'+')/M
  ##
  iSig = solve(Sig)
  lamR = eigen(R, sym=TRUE, only.val=TRUE)$val
  ## BT
  Zj = colSums(Z)
  Qb = sum(Zj*t(iSig*Zj)); lamb = R1
  pvalb = pchisq(Qb/lamb, K,lower=FALSE)
  ## MSKAT
  Q = sum(Z%*%iSig*Z); lam1 = rep(lamR, K)
  pval1 = KATpval(Q,lam1)
  ##
  L = length(rho)
  L1 = L-1; rho1 = rho[-L]
  Qw = (1-rho)*Q + rho*Qb
  pval = rep(1, L)
  pval[1] = pval1;  pval[L] = pvalb
  Lamk = vector('list', L)
  Lamk[[L]] = rep(lamb, K);  Lamk[[1]] = lam1
  tmp = sqrt(1-rho1+M*rho1) - sqrt(1-rho1)
  c1 = sqrt(1-rho1)*tmp;  c2 = tmp^2*R1/M^2
  for(k in 2:L1){
    mk = (1-rho[k])*R + c1[k]*RJ2 + c2[k]
    aak = zapsmall( abs(eigen(mk, sym=TRUE, only.val=TRUE)$val) )
    Lamk[[k]] = rep(aak[aak>0], K)
    pval[k] = KATpval(Qw[k],Lamk[[k]])
  }
  minP = min(pval)
  ## if(minP>1e-3)    return( list(p.value=c(1.5*minP,pval[c(1,L)]), pval=pval, rho.est=rho[which.min(pval)]) )
  ## p.value 
  qval = rep(0,L1)
  for(k in 1:L1) qval[k] = Liu0.qval(minP, Lamk[[k]])  ## for(k in 1:L1) qval[k] = chisum.qval(minP, Lamk[[k]])
  q1 = qchisq(minP,K,lower=FALSE)
  a1 = eigen(R, sym=TRUE);   Rh = a1$vec%*%diag(sqrt(a1$val))%*%t(a1$vec)
  Rh1 = rowSums(Rh);  H1 = outer(Rh1,Rh1)/R1
  Lamq = rep(eigen(R-R2/R1*H1, symmetric=TRUE,only.values=TRUE)$val, K)
  tauk = (1-rho1)*R2/R1 + rho1*R1
  katint = function(xpar){
    eta1 = sapply(xpar, function(eta0) min((qval-tauk*eta0)/(1-rho1)))
    KATpval(eta1,Lamq)*dchisq(xpar,K)
  }
  prec = 1e-4
  p.value = try({ minP + integrate(katint, 0,q1,  subdivisions=1e3,abs.tol=minP*prec)$val }, silent=TRUE)
  while(class(p.value)=='try-error'){
    prec = prec*2
    p.value = try({ minP + integrate(katint, 0,q1, abs.tol=minP*prec)$val }, silent=TRUE)
  }
  return(list(p.value=c(A=p.value, V=pval[1], B=pval[L]), pval=pval, rho.est=rho[which.min(pval)]) )
}



## Function to compute AT p-value as a function of the minimum p-value
MTval = function(minp, Sig, rho=0:10/10){
  K = length(rho); M = dim(Sig)[1]
  qval = rep(0,K)
  for(i in 1:K){
    if(rho[i]==0){
      qval[i] = qchisq(minp,M, lower=FALSE)
    } else if(rho[i]==1){
      qval[i] = qchisq(minp,1, lower=FALSE)
    } else{
      qval[i] = KATqval(minp, c(rep(1-rho[i],M-1),1))
    }
  }
  chiint = function(xq){
    xi = sapply(xq, function(x) min( (qval[-K]-x)/(1-rho[-K])) )
    pchisq(xi,M-1,lower=FALSE)*dchisq(xq,1)
  }
  p.val = minp + integrate(chiint, 0,qval[K])$val
  attr(p.val,'qval')  = qval
  return(p.val)
}



#' Compute the list of significant SNPs using the GWAS summary data based multi-trait association tests
#'
#' We compute the significant SNPs for three tests: minimum marginal test p-value (minP); SZ test; SZ2 test.
#'
#' @param  Z matrix of summary Z-statistics (SNPs by traits)
#' @param  Sig the estimated marginal trait correlation matrix
#' @param  alpha desired genome-wide significance level (default to 5E-8)
#' @return
#' \describe{
#'   \item{idM}{ significant SNP list for minP }
#'   \item{idS}{ significant SNP list for SZ }
#'   \item{idS2}{ significant SNP list for SZ2 }
#' }
#' @export
#' @references
#' Guo,B. and Wu,B. (2017) GWAS summary data based multi-trait SNP-set association test to identify novel genetic variants. tech rep.
Lmsatz <- function(Z, Sig, alpha=5e-8){
  if(class(Z)!='matrix') Z = as.matrix(Z)
  M = dim(Sig)[1]
  es = eigen(Sig, sym=TRUE)
  ##
  f1 = function(q0) ( 1-pmvnorm(rep(-q0,M),rep(q0,M), sigma=Sig, algorithm=Miwa(steps=256)) - alpha )*floor(1/alpha)
  q1 = qnorm(alpha/2,lower=FALSE); q2 = q1*1.1
  while( f1(q1)*f1(q2)>0 ) q2 = q2*1.1
  q0 = uniroot(f1, c(q1,q2), tol=alpha*1e-4)$root
  ##
  t1 = qchisq(alpha,1, lower=FALSE)
  t2 = KATqval(alpha,es$val,alpha*1e-4)
  ##
  SZ1 = rowSums(Z)^2/sum(Sig)
  ss1 = which(SZ1>t1)
  SZ2 = rowSums(Z^2)
  ss2 = which(SZ2>t2)
  ## Zm = abs(Z[,1]); for(i in 2:M)  Zm = pmax(Zm, abs(Z[,i])); ss3 = which(Zm>q0)
  ss3 = which(rowSums(abs(Z)>q0)>0)
  ss0 = which(rowSums(abs(Z)>q1)>0)
  ##
  tm = qchisq(alpha,M, lower=FALSE)
  Q = rowSums(Z%*%solve(Sig)*Z)
  ss4 = which(Q>tm)
  ## 
  tm2 = qchisq(alpha/2,M, lower=FALSE)
  t22 = KATqval(alpha/2,es$val,alpha*1e-4)
  ss5 = which( (Q>tm2)|(SZ2>t22) )
  ##
  return(list(idM=ss3,idS=ss1,idS2=ss2, idQ=ss4,idQ2=ss5, idm=ss0) )
}


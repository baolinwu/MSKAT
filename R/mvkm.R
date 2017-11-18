## Taken from author's website
### Required library
## library(mvtnorm)

#### Identity matrix of size dd
Identity <- function(dd){
  return(diag(1, nrow=dd))
}

###########
### Defining different kernel functions
###
### Each function taken an input as pxn matrix
###		in which each column is a sampling unit
###		and each row is one marker
###
### The output is the n-by-n kernel matrix
###########

### Linear kernel
linear.kernel <- function(X){
  return(t(X)%*%X)
}

### Quadratic kernel
quad.kernel <- function(X){
  return((1+t(X)%*%X)^2)
}

### Identical by state (IBS) kernel
IBS.kernel <- function(Z){ ## Z is the full nxm SNP matrix
  n = nrow(Z)
  K = matrix(nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in i:n) {
      K[i,j] = K[j,i] = sum(  2*(Z[j,] == Z[i,]) + (abs(Z[j,]-Z[i,])==1) )
    }
  }
  K = K/ncol(Z)
  return(K)
}



###########################
### Classical multivariate regression (parametric, Y = X\beta +  e)
###########################
Parametric_MR <- function(Y, X.list, Sigma){
### Input arguments are
###	Y :: pxn matrix, each row is one response, each column is one subject
###	X.list :: a list where each item of list is a qxn matricx (parametric design matrix),
###            where each row is one covariate (possiblly different for all items in Y)
### 	Sigma :: pxp covariance matrix (of the errors)
###
### Outout is a list containing elements
###	beta.hat :: estimate of beta (qp -by- 1 vector, first q elements correspond to first response variable and so on)
###	beta.cov :: estimated covariance matrix for estimated beta parameters
###	resid::	vector of residuals
###

  n = ncol(Y)  ## number of subjects
  p = ncol(Sigma) ## number of response variables per subject


  ## Identity of size n
  I.n = diag(1, nrow=n)

  ## big Y-vector
  Y.til = matrix(t(Y), ncol=1)

  ## inverse of sigma
  iSigma = solve(Sigma)

  ## square root of sigma-inverse
  sig.eig <- eigen(iSigma)
  isig.sqrt <- sig.eig$vectors %*% diag(sqrt(sig.eig$values)) %*% solve(sig.eig$vectors)


  ## Covariance of the entire Y vector
  Sigma.til = kronecker(Sigma, I.n)

  ## temporary var (needed later)
  iPs = kronecker(iSigma, I.n)

  ## Overall parametric design matrix (X in X-transpose-beta)
  q = nrow(X.list[[1]])
  X.mat = matrix(0, nrow=q*p, ncol=n*p)
  for (x.ind in 1:p){
    ind.set2 = c(1:n) + n*(x.ind-1)
    ind.set1 = c(1:q) + q*(x.ind-1)
    X.mat[ind.set1, ind.set2] = X.list[[x.ind]]
  }

  ## hat matrix for beta
  G = solve(X.mat %*% iPs %*% t(X.mat)) %*% X.mat %*% iPs

  ## Estimates
  beta.hat = G%*% Y.til
  beta.sd = G %*% Sigma.til %*% t(G)

  ## predictions
  pred.Y = t(X.mat)%*%beta.hat
  return(list(beta.hat = beta.hat, beta.cov = beta.sd, resid = Y.til-pred.Y))
}




#' Multivariate Kernel Machine Regression test of rare variant set
#'
#' Return the variant set association test statistic and p-value.
#'   Aadapted from Maity et. al approach for rare variant set association test using weighted linear kernel.
#' @param Y  pxn matrix, each row is one response, each column is one subject
#' @param X.list a list where each item of list is a qxn matricx (parametric design matrix),
#'            where each row is one covariate (possiblly different for all items in Y), intercept must be included!
#' @param Sigma   covariance matrix of multivariate outcomes; estimate from data if NULL
#' @param Gm    genotype matrix with SNPs in column
#' @param W    variant weight; use Beta dist weight if NULL
#' @param W.beta   Beta distribution parameters for the variant weight. Default to Beta(1,25).
#' @param n.sim   number of simulation to obtain the null distribution of the test statistic; using analytical method if NULL
#' @export
#' @references
#' Maity,A., Sullivan,P.F., Tzeng,J. (2012) Multivariate Phenotype Association Analysis by Marker-Set Kernel Machine Regression. Genet. Epidemiol. 36, 686–695.
MKMR <- function(Y, X.list, Sigma=NULL, Gm, W=NULL, W.beta=c(1,25), n.sim=NULL){
  n = ncol(Y) ## number of subjects
  p = nrow(Y)
  m = dim(Gm)[2]
  if(is.null(W)){
    W = dbeta(colMeans(Gm)/2, W.beta[1],W.beta[2])
    W = W/sum(W)*m
  } else{
    if(length(W)<=1) W = rep(1,m)
  }
  Gm = t(t(Gm)*W)

  ## Identity of size n
  I.n = diag(1, nrow=n)

  ## big Y-vector
  Y.til = matrix(t(Y), ncol=1)

  ## inverse of sigma
  if(is.null(Sigma)){
    res = matrix(0, n,p)
    for(k in 1:p){
      res[,k] = lm(Y[k,] ~ t(X.list[[k]]))$res
    }
    Sigma = cov(res)
  }
  iSigma = solve(Sigma)

  ## square root of sigma-inverse
  sig.eig <- eigen(iSigma)
  isig.sqrt <- sig.eig$vectors %*% diag(sqrt(sig.eig$values)) %*% solve(sig.eig$vectors)


  ## Covariance of the entire Y vector
  iSigma.til = kronecker(iSigma, I.n)

  ## overall Kernel matrix
  K.mat = matrix(0, ncol=n*p, nrow=n*p)
  Kg = Gm%*%t(Gm)
  for (k.ind in 1:p) {
    ind.set = 1:n + n*(k.ind-1)
    K.mat[ind.set, ind.set] = Kg ## K.list[[k.ind]]
  }

  ## Overall parametric design matrix (X in X-transpose-beta)
  q = nrow(X.list[[1]])
  X.mat = matrix(0, nrow=q*p, ncol=n*p)
  for (x.ind in 1:p) {
    ind.set2 = c(1:n) + n*(x.ind-1)
    ind.set1 = c(1:q) + q*(x.ind-1)
    X.mat[ind.set1, ind.set2] = X.list[[x.ind]]
  }

  ## Hat matrix for beta
  G = solve(X.mat %*% iSigma.til %*% t(X.mat)) %*% X.mat %*% iSigma.til

  ## estimate under H0
  beta.hat = G%*% Y.til

  ## prediction of Y (under null)
  pred.Y = t(X.mat)%*%beta.hat

  ## Test statistics
  test.stat.our = t(Y.til-pred.Y) %*%    iSigma.til%*%K.mat%*%iSigma.til      %*% (Y.til-pred.Y)

  ## Simulation based test and p-value
  P0 = iSigma.til - iSigma.til%*%t(X.mat)%*%G
  P0.eig = eigen(P0)
  vtmp = P0.eig$values * (P0.eig$values > 0)
  P0.sqrt = P0.eig$vectors %*% diag(sqrt(vtmp)) %*% solve(P0.eig$vectors)

  M.mat = P0.sqrt %*% K.mat %*% P0.sqrt
  M.eig = matrix(Re(eigen(M.mat)$values), ncol=1)

  if(!is.null(n.sim)){
    sim.tmp = rep(NA, n.sim)

    sim.r = matrix(rchisq(n.sim*nrow(M.eig), df=1), ncol=nrow(M.eig))
    sim.tmp = sim.r %*% M.eig

    pv = mean(c(sim.tmp) > c(test.stat.our))
  } else{
    pv = KAT.pval(test.stat.our, M.eig)
  }
  return(list(test.stat = c(test.stat.our), p.value = pv))
}


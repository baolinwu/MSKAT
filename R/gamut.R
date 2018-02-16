#' GAMuT: Gene Association with Multiple Traits 
#'
#' Return the multi-trait variant set association test p-values from the Broadaway et. al (2016) approach.
#'  Both tests from linear kernel and projection matrix are computed.
#' @param PHENO  outcome matrix. sample in rows.
#' @param GENO   genotype matrix. sample in rows.
#' @param COVAR  covariate matrix to be adjusted. sample in rows.
#' @param W    variant weight; use Beta dist weight if NULL
#' @param W.beta   Beta distribution parameters for the variant weight. Default to Beta(1,25).
#' @return a p-value vector p.value with two components
#' \describe{
#'   \item{proj}{ p-value for GAMuT with the projection matrix kernel }
#'   \item{linear}{ p-value for GAMuT with the linear kernel }
#' }
#' @export
#' @references
#' Broadaway, K.A., Cutler, D.J., Duncan, R., Moore, J.L., Ware, E.B., Jhun, M.A., Bielak, L.F., Zhao, W., Smith, J.A., Peyser, P.A., Kardia, S.L.R., Ghosh, D., Epstein, M.P. (2016). A Statistical Approach for Testing Cross-Phenotype Effects of Rare Variants. \emph{The American Journal of Human Genetics}, 98, 525--540.
#'
#' Original R codes posted at \url{http://genetics.emory.edu/labs/epstein/software/gamut}
gamut <- function(PHENO, GENO, COVAR=NULL, W=NULL, W.beta=c(1,25)){
  ##===============================================================================
  ## STEP 2:  residualize the phenotypes in PHENO on the covariates in COVAR
  ## Note: Each phenotype in PHENO will be residualized on each covariate in COVAR
  ##===============================================================================
  ##------------------------------------------------
  ## If you wish to only residualize phenotypes on a subset of covariates in COVAR, 
  ## make sure to make appropriate changes to the code below
  ##------------------------------------------------
  if(is.null(COVAR)){
    residuals = PHENO ## add the option of no covariates
  } else{
    residuals = matrix(NA, nrow=nrow(PHENO), ncol=ncol(PHENO)) # Matrix of residualized phenotypes
    for(l in 1:ncol(PHENO)){
      model.residual = lm(PHENO[,l] ~ as.matrix(COVAR))
      residuals[,l] = resid(model.residual)
    }
  }
  
  ##===============================================================================
  ## STEP 3:  form the phenotypic similarity matrix Yc (using notation from the AJHG paper)
  ## and corresponding eigenvalues lambda_Y based on residualized phenotypes from Step 2
  ##===============================================================================
  ## centered & scaled matrix of residualized phenotypes:
  P0 <- apply(residuals, 2, scale, scale=T, center=T)
  P  <- as.matrix(P0) # convert dataframes into matrix
  
  ##------------------------------------------------
  ## construct Yc using the projection matrix
  ##------------------------------------------------
  ## function for constructing the projection matrix and corresponding eigenvalues:
  proj_pheno = proj_GAMuT_pheno(P)
  Yc1 = proj_pheno$Kc                # Projection matrix
  lambda_Y1 = proj_pheno$ev_Kc       # Eigenvalues of Yc
  
  ## function for constructing the linear phenotype kernel and eigenvalues:
  linear_pheno <- linear_GAMuT_pheno(P) 
  Yc2 = linear_pheno$Kc           # Linear kernel similarity matrix
  lambda_Y2 = linear_pheno$ev_Kc  # Eigenvalues of Yc


  ##===============================================================================
  ## STEP 4:  form the genotypic similarity matrix Xc (using notation from the AJHG paper)
  ## and corresponding eigenvalues lambda_X
  ##===============================================================================
  
  ## We assume weighted linear kernel for genotypes where weights are function of 
  ## minor-allele frequency (MAF). We apply the beta-distribution MAF weights of 
  ## Wu et al. (2011) in the SKAT paper
  
  ## form the variant weights:
  MAF = colMeans(GENO)/2                      # sample MAF of each variant in the sample
  if(is.null(W)){
    beta_weight = dbeta(MAF,W.beta[1],W.beta[2]) ## /dbeta(0,W.beta[1],W.beta[2]) # assume beta-distribution weights
    beta_weight = beta_weight/sum(beta_weight)*ncol(GENO) ## avoid numerical errors
    ## Note: one can use other weight functions by simply recoding beta_weight as desired
  } else{
    beta_weight = W/sum(W)*ncol(GENO)
  }
  ## Weighted rare variants
  G0 = as.matrix(GENO)%*%diag(beta_weight) # Weighted rare variants # corrected 1
  ## G0 = beta_weight*GENO ### original 1 ## wrong
  ## G0 = t(beta_weight*t(GENO)) ## corrected
  G  = as.matrix(scale(G0,center=T,scale=F)) # Centered genotype matrix

  ## function for constructing the weighted linear kernel matrix for genotypes and eigenvalues:
  linear_geno <- linear_GAMuT_geno(G) 
  Xc <- linear_geno$Lc            # Linear kernel similarity matrix
  lambda_X <- linear_geno$ev_Lc   # Eigenvalues of Xc
  ##===============================================================================
  ## STEP 5: Construct the GAMuT Test and obtain the p-value
  ##===============================================================================
  ## pval = rep(NA,2); names(pval) = c('proj', 'linear')
  pval1 = TestGAMuT(Yc1,lambda_Y1,Xc,lambda_X)
  pval2 = TestGAMuT(Yc2,lambda_Y2,Xc,lambda_X)
  return( p.value=c('proj'=pval1, 'linear'=pval2) )
}


## GAMuT-functions.R
##
## 2016-March-13
##
## this script contains functions used in the main program
## for GAMuT analysis
##
##----------------------------------------------------------------------
## descriptions of individual functions:
##----------------------------------------------------------------------
##
## * proj_GAMuT_pheno
##   constructs projection matrix and corresponding eigenvalues for phenotypes
##
## * linear_GAMuT_pheno
##   constructs linear kernel matrix and corresponding eigenvalues for phenotypes
##
## * linear_GAMuT_geno
##   constructs linear kernel matrix and corresponding eigenvalues for genotypes 
##
## * TestGAMuT
##   constructs GAMuT statistic and returns p-value


##----------------------------------------------------------------------
## phenotypic similarity:  projection matrix
##----------------------------------------------------------------------
proj_GAMuT_pheno <- function(X){
    out = vector("list", 2)
    names(out) = c("Kc", "ev_Kc")

    out$Kc = X %*% solve(t(X) %*% X) %*% t(X)   # projection matrix
    out$ev_Kc = rep(1, ncol(X))                 # find eigenvalues
    return(out)	
}

##----------------------------------------------------------------------
## phenotypic similarity:  linear kernel
##----------------------------------------------------------------------
linear_GAMuT_pheno <- function(X){
    out = vector("list", 2)
    names(out) = c("Kc", "ev_Kc")
    
    Kc = t(X) %*% X   # transposed kernel to find eigenvalues
    ev_Kc = eigen(Kc, symmetric=T, only.values=T)$values  
    
    out$Kc = X %*% t(X) # similiarity kernel for test
    out$ev_Kc = ev_Kc[ev_Kc > 1e-08]
    return(out)
}

##----------------------------------------------------------------------
## genotypic similarity:  linear kernel
##----------------------------------------------------------------------
linear_GAMuT_geno <- function(X){
    out = vector("list", 2)
    names(out) = c("Lc", "ev_Lc")
    
    Lc = t(X) %*% X   # transposed kernel to find eigenvalues
    ev_Lc = eigen(Lc, symmetric=T, only.values=T)$values  
    
    out$Lc = X %*% t(X) # similiarity kernel for test
    out$ev_Lc = ev_Lc[ev_Lc > 1e-08]
    return(out)	 
}


##----------------------------------------------------------------------
## constructing the GAMuT statistic and deriving p-value:
##----------------------------------------------------------------------
TestGAMuT <- function(Yc, lambda_Y, Xc, lambda_X) {

    ## test statistic:
    m = nrow(Yc) # number of subjects in study
    GAMuT = (1/m) * sum(sum(t(Yc) * Xc))  
    

    ## populate vector of all pairwise combination of eigenvalues
    ## from the phenotype and genotype similarity matrices:
    Z <- (as.matrix(lambda_Y)) %*% t(as.matrix(lambda_X))
    Zsort <- sort(Z, decreasing=T)

    ## derive p-value of GAMuT statistic:
  scoredavies = GAMuT*m^2
  results_score <- davies(scoredavies, Zsort)
  davies_pvalue <- (results_score$Qq)

   return(davies_pvalue)
  ## davies_pvalue1 = KAT.pval(scoredavies, Zsort)
  ##  return(c(davies_pvalue, davies_pvalue1))
} 

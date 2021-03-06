% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PheSATS.R
\name{PheSATS}
\alias{PheSATS}
\title{GWAS summary data based phenome-wide SNP-set association test}
\usage{
PheSATS(Z, Sig, R)
}
\arguments{
\item{Z}{M by K matrix of summary Z-statistics for M variants across K traits}

\item{Sig}{estimated trait correlation matrix (K by K)}

\item{R}{variant LD correlation matrix (M by M)}
}
\value{
p-value vector: FCT, SBT, S2T, TCT, and S2E (residual effects excluding SBT)
}
\description{
Efficient and powerful summary data based phenome-wide SNP-set association tests: quadratic test (S2T), burden test (SBT), their adaptive test:
Fisher combination test (FCT), and Tippett combination test (TCT).
}
\examples{
K = 4; M = 20
R = cor(matrix(rnorm(500*M),500,M)*sqrt(0.8)+rnorm(500)*sqrt(0.2))
Y = matrix(rnorm(100*K), 100,K)*sqrt(0.8)+rnorm(100)*sqrt(0.2); Y[,1] = -Y[,1]
Sig = cor(Y)
Rc = t(chol(R)); Sc = chol(Sig)
Z0 = Rc\%*\%matrix(rnorm(K*M),M,K)\%*\%Sc
PheSATS(Z0,Sig,R)
js = sample(1:M, size=round(M*0.4))
Z1 = Z0; Z1[js,] = Z1[js,] + 1.5
PheSATS(Z1,Sig,R)
ij = cbind(sample(1:M, size=20, rep=TRUE), sample(1:K, size=20, rep=TRUE))
Z2 = Z0; Z2[ij] = Z2[ij] + rnorm(20)*1.5
PheSATS(Z2,Sig,R)
}
\references{
Guo,B., Bishop,J., Liu,N. and Wu,B. (2018) Phenome-wide SNP-set association test based on GWAS summary data to identify novel disease-gene association. tech rep.
}

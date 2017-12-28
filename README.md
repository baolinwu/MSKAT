# MSKAT
 - An R pacakge implementing various statistical methods for testing multi-trait variant-set association

-----
## Testing the rare variant set association with multiple continuous phenotypes
 - Two multivariate versions of SKAT are implemented.
 - A variation of the multivariate kernel machine regression (MKMR) is adapted to rare variant set association test.
    - Ref: Maity,A., Sullivan,P.F., Tzeng,J. (2012) Multivariate Phenotype Association Analysis by Marker-Set Kernel Machine Regression. *Genet. Epidemiol.*, 36, 686–695.
 - Multi-trait SKAT
    - Ref: Wu,B., Pankow,J.S. (2016). Sequence kernel association test of multiple continuous phenotypes. *Genetic Epidemiology*, 40(2), 91-100.

 - Sample R codes
```r
library(MSKAT)
Y = matrix(rnorm(2000), 1000,2)
X = matrix(rnorm(2000), 1000,2)
G = matrix(rbinom(10000,2,0.02), 1000,10)
MSKAT(MSKAT.cnull(Y,X), G, W.beta=c(1,25))
## Time consuming
X.list = vector('list',2)
for(i in 1:2) X.list[[i]] = rbind(1, t(X))
MKMR(t(Y),X.list, Gm=G, W.beta=c(1,25))
```



-----
## Multi-trait SNP-set association tests using GWAS summary data
 - Reference
    - Guo,B. and Wu,B. (2017) Powerful and efficient SNP-set association tests across multiple phenotypes using GWAS summary data. tech rep. 
 - Efficient and power MSATS tests: variance components test (VC), burden type test (BT), adaptive test (AT)
    - All test p-values are efficiently and analytically computed: no need of Monte Carlo sampling and extremely scalable. 
    - AT has robust performance and truly combines the strength of both VC and BT.
 - Sample R codes
```r
K = 4; M = 20
R = cor(matrix(rnorm(500*M),500,M)*sqrt(0.8)+rnorm(500)*sqrt(0.2))
Y = matrix(rnorm(100*K), 100,K)*sqrt(0.8)+rnorm(100)*sqrt(0.2)
Y[,1] = -Y[,1]; Sig = cor(Y)
Z = matrix(rnorm(K*M),M,K)
js = sample(1:M, size=round(M*0.4))
Z[js,] = Z[js,] + 1.5
MSATS(Z,Sig,R)
Z = matrix(rnorm(K*M),M,K)
ij = cbind(sample(1:M, size=20, rep=TRUE), sample(1:K, size=20, rep=TRUE))
Z[ij] = Z[ij] + rnorm(20)*1.5
MSATS(Z,Sig,R)
```




-----
## Testing the rare variant set association across multiple traits.

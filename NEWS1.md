# MSKAT
 - An R package implementing various statistical methods for testing multi-trait variant-set association

-----
## SNP-set association tests across multiple phenotypes using GWAS summary data
 - Multi-trait SNP-set Association Tests using GWAS Summary data (MSATS)
 - Efficient and powerful MSATS tests: variance components test (VC), burden type test (BT), adaptive test (AT)
    - All test p-values are efficiently and analytically computed: extremely scalable without the need of Monte Carlo sampling.
    - The AT has robust performance and can truly combine the strength of both VC and BT.
 - Reference
    - Guo,B. and Wu,B. (2017) Powerful and efficient SNP-set association tests across multiple phenotypes using GWAS summary data. tech rep. 
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

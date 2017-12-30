# MSKAT
 - An R package implementing various statistical methods for testing multi-trait variant-set association

-----
## Testing the rare variant set association with multiple continuous phenotypes
 - Two multivariate versions of SKAT (MSKAT) are implemented.
 - Reference
    - Wu,B. and Pankow,J.S. (2016). Sequence kernel association test of multiple continuous phenotypes. *Genetic Epidemiology*, 40(2), 91-100.
 - A variation of the multivariate kernel machine regression (MKMR) is also adapted to rare variant set association test.
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
## Testing the rare variant set association across multiple traits.

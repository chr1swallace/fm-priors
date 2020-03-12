### use real medium LD region to simulate a single CV fine mapping region

setwd("~/Projects/fm-priors")
source("common.R")
library(simGWAS)
library(rmutil) # laplace distribution
library(magrittr)
library(data.table)
library(randomFunctions)
library(corrcoverage)

args <- getArgs(defaults=list(NN=5000,# number of cases = number of controls
                              NREP=5, # replicates per scenario
                              cvtype="friendly"), # or lonely
                numeric=c("NN","NREP"))

tmp <- readRDS("medium.RDS") # haplotype frequencies and LD matrix

h <- tmp$h
maf <- colMeans(h)
LD <- tmp$LD

freq <- as.data.frame(h+1)
freq$Probability <- 1/nrow(freq)
snps <- colnames(freq)[-ncol(freq)]

## analysis functions

calc_bf <- function(W,bhat,V) {
  h0 <- dnorm(bhat,sd=sqrt(V),log=TRUE)
  n <- dnorm(bhat,sd=sqrt(V+W),log=TRUE)
  lambda <- W2l(W)
  f <- function(x,b,v)
    dnorm(b, mean=x, sd=sqrt(v)) * dlaplace(x, 0, lambda)
  fint <- function(b,v)
    integrate(f,lower=0,upper=Inf,b=b,v=v)$value
  vfint <- Vectorize(fint)
  l <- matrix(vfint(bhat,V),nrow(bhat),ncol(bhat))
  list(n=n-h0, l=l-h0)
}

calc_pp <- function(bf) { # rows=replicates, cols=snps
  rs <- log(rowSums(exp(bf)))
  exp(bf - matrix(rs, nrow(bf), ncol(bf)))
}

summ <- function(pp,iCV) {
  r <- apply(pp, 1, rank)  %>% matrix(., nrow(pp), ncol(pp))
  s=data.frame(PP=pp[,iCV], rank=r[,iCV])
}

simmer <- function(iCV,truedn=c("norm","laplace"),trueW=0.2^2) {
  dn <- match.arg(truedn)
  ## sample beta from prior
  beta <- switch(dn,
                 "norm"=rnorm(1, 0, sqrt(trueW)),
                 "laplace"=rlaplace(1,0,W2l(trueW)))
  zsim = simulated_z_score(N0=args$NN, # number of controls
                           N1=args$NN, # number of cases
                           snps=snps,
                           W=snps[iCV], # causal variants, subset of snps
                           gamma.W=beta, # log odds ratios
                           freq=freq,nrep=args$NREP)
  vbetasim <- simulated_vbeta(N0=10000, # number of controls
                              N1=10000, # number of cases
                              snps=snps, # column names in freq of SNPs for which Z scores should be generated
                              W=snps[iCV], # causal variants, subset of snps
                              gamma.W=beta, # log odds ratios
                              freq=freq, # reference haplotypes
                              nrep=args$NREP)
  betasim <- zsim * sqrt(vbetasim)
  ## Wval <- c(W1=0.1,W2=0.2,W3=0.3)
  BF <- calc_bf(W=0.2^2,bhat=betasim, V=vbetasim) # always analyse with default W
  PP <- lapply(BF, calc_pp)
  PPsumm <- lapply(PP, summ, iCV=iCV)
  res <- cbind(as.data.table(PPsumm), beta=beta, truedn=dn, trueW=trueW, z=zsim[,iCV])
}


## scenarios
## sample a CV with lots of LD friends
cases <- expand.grid(trueW=Wvec,truedn=c("norm","laplace"),stringsAsFactors=FALSE)
nfriends <- apply(LD^2>0.5,1,sum)
iCV <- switch(args$cvtype,
              friendly=sample(which(nfriends>30 & pmin(maf,1-maf)>0.05), 1),
              lonely=sample(which(nfriends==1 & pmin(maf,1-maf)>0.05), 1))

for(j in 1:200) {
  res <- lapply(1:nrow(cases), function(i) {
    cat(i,"\t")
    simmer(iCV,truedn=cases$truedn[i],trueW=cases$trueW[i])
  })  %>% rbindlist()
  res$cvtype <- factor(args$cvtype,levels=c("friendly","lonely"))
  res$nn <- args$NN
  
  f <- tempfile(tmpdir="~/scratch/fm-priors",fileext=".csv.gz")
  fwrite(res,f)
}



setwd("~/Projects/fm-priors")
source("common.R")
library(rmutil) # laplace distribution
library(ggplot2)
library(data.table)
library(magrittr)
library(hrbrthemes); theme_set(theme_ipsum())

## analytical marginal likelihood
marg <- function(W,data) {
  n0 <- dnorm(data$bhat,sd=sqrt(data$V),log=TRUE)
  n <- dnorm(data$bhat,sd=sqrt(data$V+W),log=TRUE)
  lambda <- W2l(W)
  f <- function(x)
    dnorm(data$bhat, mean=x, sd=sqrt(data$V)) * dlaplace(x, 0, lambda)
  l2 <- integrate(f, -Inf, Inf)
  c(unlist(data), W=W,h0=n0,h1.n=n,h1.l=log(l2$value))
}

simmer <- function(truedn=c("norm","laplace"),trueW=0.2^2) {
  dn <- match.arg(truedn)
  ## sample beta from prior
  beta <- switch(dn,
                 "norm"=rnorm(1, 0, sqrt(trueW)),
                 "laplace"=rlaplace(1,0,W2l(trueW)))
  ## sample bhat | beta
  V <- sample(c(0.001,0.01,0.1),1)
  bhat <- rnorm(1, beta, sqrt(V))
  ## analyse
  data <- list(V=V,bhat=bhat)
  ret <- lapply(Wvec, marg, data)  %>% do.call("rbind",.)  %>% as.data.frame()
  ret$trueW <- trueW
  ret$truedn <- dn
  ret$beta <- beta
  ret
}

simmer()
n <- replicate(1000,simmer(),simplify=FALSE)  %>% rbindlist()
l <- replicate(1000,simmer("laplace"),simplify=FALSE)  %>% rbindlist()
results <- rbind(n,l)
results[,bf.n:=h1.n - h0]
results[,bf.l:=h1.l - h0]
results[,z:=abs(bhat/sqrt(V))]

## considerable differences in marginal likelihood
## rather less on Bayes factor scale
ymn <- min(c(results$h1.n,results$h1.l))
ymx <- max(c(results$h1.n,results$h1.l))
ymn2 <- min(c(results$bf.n,results$bf.l))
ymx2 <- max(c(results$bf.n,results$bf.l))
library(cowplot)
p1 <- ggplot(results,aes(x=z,y=h0)) + geom_point() + ggtitle("Marg. lhood","H0") + ylab("log likelihood")
p2 <- ggplot(results,aes(x=z,y=h1.n,col=factor(W),pch=factor(V))) + geom_point() + ggtitle("Marg. lhood","HA, normal") +
  scale_colour_ipsum() +
  scale_shape(name="Simulation V",solid=FALSE) +
  ## scale_colour_manual("Analysis W",values=tableau_colours[1:3]) +
  ylim(ymn,ymx) + ylab("log likelihood")
p3 <- ggplot(results,aes(x=z,y=h1.l,col=factor(W),pch=factor(V))) + geom_point() + ggtitle("Marg. lhood","HA, Laplace") +
  scale_colour_ipsum("Analysis W") +
  scale_shape(name="Simulation V",solid=FALSE) +
  ylim(ymn,ymx) + ylab("log likelihood")
p4 <- ggplot(results,aes(x=z,y=h1.n-h0,col=factor(W))) + geom_point() + ggtitle("BF", "normal") + geom_smooth(se=FALSE) +
  scale_colour_ipsum() +
  ## scale_shape(name="Simulation V",solid=FALSE) +
  ylim(ymn2,ymx2) + ylab("log Bayes factor")
p5 <- ggplot(results,aes(x=z,y=h1.l-h0,col=factor(W))) + geom_point() + ggtitle("BF", "Laplace") + geom_smooth(se=FALSE) +
  scale_colour_ipsum("W") +
  ## scale_shape(name="Simulation V",solid=FALSE) +
  ylim(ymn2,ymx2) + ylab("log Bayes factor")
# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  p3 + theme(legend.box.margin = margin(0, 0, 0, 12))
)
plots <- list(p1,p2,p3,p4,p5)  %>% lapply(., function(p) p + theme(legend.position="none"))
plots <- c(plots,list(legend))
plot_grid(plotlist=plots)
ggsave("marg-lhood-bf.png",height=8,width=8)

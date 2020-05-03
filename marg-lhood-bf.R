setwd("~/Projects/fm-priors")
source("common.R")
library(rmutil) # laplace distribution
library(ggplot2)
library(data.table)
library(magrittr)
## library(hrbrthemes); theme_set(theme_ipsum())
## remotes::install_github("chr1swallace/seaborn_palettes")
library(seaborn)
set.seed(42) # reproducibility

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
## rather smaller differences in Bayes factor
ymn <- min(c(results$h1.n,results$h1.l))
ymx <- max(c(results$h1.n,results$h1.l))
ymn2 <- min(c(results$bf.n,results$bf.l))
ymx2 <- max(c(results$bf.n,results$bf.l))
library(cowplot)
## marg under H0
p1 <- ggplot(results,aes(x=z,y=h0)) +
  geom_point() +
  ggtitle("","H0") +
  ylab("log likelihood") +
  theme_cowplot(12) + 
  background_grid(major = "xy", minor = "none") +
  theme(plot.title = element_text(color="black", size=14, face="plain"))

## marg under HA, normal
p2 <- ggplot(results,aes(x=z,y=h1.n,col=factor(W),pch=factor(V))) + geom_point() + ggtitle("","HA, Gaussian") +
  scale_shape(name="Simulation V",solid=TRUE) +
  scale_colour_seaborn("Analysis W") +
  ## scale_shape(name="Simulation V",solid=FALSE) +
  ## scale_colour_manual("Analysis W",values=tableau_colours[1:3]) +
  ylim(ymn,ymx) + ylab("log likelihood") +
  theme_cowplot(12) + 
  background_grid(major = "xy", minor = "none") +
  theme(plot.title = element_text(color="black", size=14, face="plain")) 

p3 <- ggplot(results,aes(x=z,y=h1.l,col=factor(W),pch=factor(V))) +
  geom_point() +
  ggtitle("","HA, Laplace") +
  labs(col = "Analysis W") +
  scale_shape(name="Simulation V",solid=TRUE) +
  scale_colour_seaborn("Analysis W") +
  ylim(ymn,ymx) + ylab("log likelihood") +
  theme_cowplot(12) + 
  background_grid(major = "xy", minor = "none")+
  theme(plot.title = element_text(color="black", size=14, face="plain"))

p4 <- ggplot(results,aes(x=z,y=h1.n-h0,col=factor(W))) +
  ggtitle("", "Gaussian") +
  geom_smooth(se=FALSE) +
  geom_point() +
  scale_colour_seaborn("Analysis W") +
  ylim(ymn2,ymx2) + ylab("log BF")+
  theme_cowplot(12) + 
  background_grid(major = "xy", minor = "none")+
  theme(plot.title = element_text(color="black", size=14, face="plain"))

p5 <- ggplot(results,aes(x=z,y=h1.l-h0,col=factor(W))) +
  geom_smooth(se=FALSE) +
  geom_point() +
  ggtitle("", "Laplace") +
  scale_colour_seaborn("Analysis W") +
  ylim(ymn2,ymx2) + ylab("log BF")+
  theme_cowplot(12) + 
  background_grid(major = "xy", minor = "none")+
  theme(plot.title = element_text(color="black", size=14, face="plain"))

## extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  p3 + theme(legend.box="horizontal",
             legend.box.background = element_rect(colour = "black"),
             legend.box.margin = margin(6, 6, 6, 6)) #legend.box.margin = margin(0, 0, 0, 12))
)

plots <- list(p1,p2,p3,p4,p5)  %>% lapply(., function(p) p + theme(legend.position="none"))
plots <- c(plots,list(legend))

plus <- ggdraw() +
  draw_label("+", x = 0.5, y=0.5,size=15, hjust = 0.5)
minus <- ggdraw() +
  draw_label("-", x = 0.5, y=0.55,size=30, hjust = 0.5)
equals <- ggdraw() +
  draw_label("=", x = 0.5, y=0.5,size=30, hjust = 0.5)
bf <- ggdraw() +
  draw_label("log BF", x = 0.5, y=0.5,size=15, hjust = 0.5)
H1 <- ggdraw() +
  draw_label("log marginal likelihood, HA", x = 0.5, y=0.5,size=15, hjust = 0.5)
H0 <- ggdraw() +
  draw_label("log marginal likelihood, H0", x = 0.5, y=0.5,size=15, hjust = 0.5)

w <- 0.1
toprow <- plot_grid(H1,minus,H0,equals,bf,nrow=1,rel_widths=c(1,w,1,w,1))+
  theme(plot.background = element_rect(color = "black", size = 3))
botright <- plot_grid(plots[[4]],plots[[5]],nrow=2)
botleft <- plot_grid(plots[[2]],plots[[3]],nrow=2)
botbotmid <- plot_grid(NULL,legend,nrow = 1, rel_widths = c(0.8,2))
botmid <- plot_grid(NULL,plots[[1]],botbotmid,nrow=3,rel_heights=c(1,2,1))
bot <- plot_grid(botleft,NULL,botmid,NULL,botright,nrow=1,rel_widths=c(1,w,1,w,1))

final <- plot_grid(toprow,bot,nrow=2,rel_heights=c(0.1,0.9))

# ggsave("marg-lhood-bf.png",plot=final,height=8,width=10)
ggsave("marg-lhood-bf.png",plot=final,height=8,width=8)

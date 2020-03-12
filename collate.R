setwd("~/Projects/fm-priors")
library(simGWAS)
library(rmutil) # laplace distribution
library(magrittr)
library(data.table)
library(randomFunctions)
library(ggplot2)
library(hrbrthemes); theme_set(theme_ipsum())

res <- fread(cmd="zcat ~/scratch/fm-priors/*.csv.gz|grep -F -v PP|grep -F -v NA")
f <- list.files("~/scratch/fm-priors",full=TRUE)[1]
nm <- scan(f,what="",nlines=1,sep=",")
setnames(res,nm)

head(res)
nrow(res)
res <- res[abs(z)<21]
table(res$truedn, res$trueW)
table(res$cvtype, res$nn)

################################################################################
## dn of PP/rank
m <- melt(res[],c("cvtype","nn","beta","truedn","trueW","z"),
          list(pp=grep("PP",names(res),value=TRUE),
               rank=grep("rank",names(res),value=TRUE)))
vlev <- grep("PP",names(res),value=TRUE)  %>% gsub("W|.PP","",.)
vlev
m[,dn:=vlev[variable]]

head(m)
m$truedn <- paste0("Simulation: ",m$truedn)
## m$trueW <- paste0("Simulation: W=",m$trueW)

m[,cond:=paste(cvtype,nn,trueW)]
m[,nn:=factor(nn, levels=c(2000,5000,10000))]
save(m, file="pp-sims.RData")

## can restart from here if needed

library(cowplot)

## key
k <- melt(unique(m[,.(trueW,cond,cvtype,nn)]),c("cond"))
k[,c("cvtype","nn","trueW"):=tstrsplit(cond, " ")]
k[,nn:=factor(nn, levels=c(2000,5000,10000))]
k$value <- factor(as.character(k$value) %>%
                  sub("friendly","high",.) %>%
                  sub("lonely","low",.),
                  levels=c("2000","5000","10000","low","high",Wvec))
k$variable %<>% sub("cvtype","LD",.)  %>% sub("nn","N",.)  %>% sub("trueW","W0",.)
k <- k[order(value)]
head(k)

levcond <- outer(c(2000,5000,10000),c("friendly","lonely"),function(x,y) paste(y,x))  %>% as.vector()  %>%
  outer(Wvec,., function(x,y) paste(y,x))  %>% as.vector()
levcond

k$cond <- factor(as.character(k$cond),levels=levcond)
m$cond <- factor(as.character(m$cond),levels=levcond)

p1 <- ggplot(k, aes(x=cond,y=value,col=variable)) +
  ## geom_hline(aes(yintercept=value)) +
  geom_point(size=3) +
  facet_grid(variable ~ cvtype + nn, scales="free", switch="y") +
  scale_colour_manual(values=tableau_colours[1:3]) +
  theme(axis.text.x=element_blank(),
        legend.position="none",
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        strip.text.y=element_text(hjust=1,family=""),
        strip.text.x=element_blank(),
        strip.background=element_blank(),
        strip.placement="outside")
p1

m[dn=="l",dn:="Laplace"]
m[dn=="n",dn:="Gaussian"]
m[,truedn:=sub("laplace","Laplace",truedn)]
m[,truedn:=sub("norm","Gaussian",truedn)]
m[,x:=as.numeric(cond) + ifelse(dn=="Gaussian",-0.2,0.2)]
p0 <- ggplot(m, aes(x=cond,y=pp,col=dn)) +
  ## geom_boxplot(outlier.shape=NA,notch=TRUE)  +
  ## geom_split_violin() +
##   stat_summary(position = position_dodge(width = .5),
##                fill="white",
##                geom="crossbar", width=0.15) +
  stat_summary(fun.ymin = function(z) { quantile(z,0.01) },
               fun.ymax = function(z) { quantile(z,0.99) },
               position = position_dodge(width = .5),
               fill="white",
               ## fun.y = median,
               fun.y = mean,
## fun.data="median_iqr", fun.args = list(mult=1), 
               geom="linerange") +
   stat_summary(fun.ymin = function(z) { quantile(z,0.1) },
               fun.ymax = function(z) { quantile(z,0.9) },
               position = position_dodge(width = .5),
               fill="white",
               ## fun.y = median,
               fun.y = mean,
## fun.data="median_iqr", fun.args = list(mult=1), 
               geom="crossbar", width=0.20) +
  facet_grid( truedn ~ cvtype + nn, switch="y",scales="free_x") +
  scale_colour_manual("Analysis distribution",values=tableau_colours[c(8:9)]) +
  scale_fill_manual("Analysis distribution",values=tableau_colours[c(8:9)]) +
  ## scale_colour_ipsum() +
  ## scale_fill_ipsum() +
  ## theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0))
  theme(axis.text.x=element_blank(),
        legend.position="top",
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        strip.text.y=element_text(hjust=1,family=""),
        strip.text.x=element_blank(),
        strip.background=element_blank(),
        strip.placement="outside")
## p0
plot_grid(p0,p1,ncol=1,rel_heights=c(0.7,0.3), align="v",axis="b") 

ggsave("pp-sims.png",height=8,width=8)


par(mfrow=c(1,3))
f <- function(W) 
  rnorm(10000,0,sd=sqrt(W))  %>% abs()  %>% exp()  %>% quantile(., c(0.5, 0.75,0.9, 0.99))

lapply(c(0.1,0.15,0.2,0.25,0.3)^2, f)
lapply(c(0.01,0.02,0.03,0.04,0.05),f)

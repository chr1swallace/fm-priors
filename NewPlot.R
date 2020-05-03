library(cowplot)

p1 <- ggplot(results,aes(x=z,y=h0)) + geom_point(cex = 0.5) + ggtitle("H0") + 
  ylab("log likelihood") + 
  theme_cowplot(12) + 
  background_grid(major = "xy", minor = "none") +
  theme(plot.title = element_text(color="black", size=14, face="plain"))

p2 <- ggplot(results,aes(x=z,y=h1.n,col=factor(W),pch=factor(V))) + geom_point(cex = 2) + ggtitle("HA, Gaussian") +
  # scale_colour_ipsum() +
  scale_shape(name="Simulation V",solid=TRUE) +
  ## scale_colour_manual("Analysis W",values=tableau_colours[1:3]) +
  ylim(ymn,ymx) + ylab("log likelihood") +
  theme_cowplot(12) + 
  background_grid(major = "xy", minor = "none") +
  theme(plot.title = element_text(color="black", size=14, face="plain"))


p3 <- ggplot(results,aes(x=z,y=h1.l,col=factor(W),pch=factor(V))) + geom_point(cex = 2) + ggtitle("HA, Laplace") +
  labs(col = "Analysis W") +
  scale_shape(name="Simulation V",solid=TRUE) +
  ylim(ymn,ymx) + ylab("log likelihood") +
  theme_cowplot(12) + 
  background_grid(major = "xy", minor = "none")+
  theme(plot.title = element_text(color="black", size=14, face="plain"))

p4 <- ggplot(results,aes(x=z,y=h1.n-h0,col=factor(W))) + geom_point() + ggtitle("Gaussian") + geom_smooth(se=FALSE) +
  # scale_colour_ipsum() +
  ## scale_shape(name="Simulation V",solid=FALSE) +
  ylim(ymn2,ymx2) + ylab("log BF")+
  theme_cowplot(12) + 
  background_grid(major = "xy", minor = "none")+
  theme(plot.title = element_text(color="black", size=14, face="plain"))


p5 <- ggplot(results,aes(x=z,y=h1.l-h0,col=factor(W))) + geom_point() + ggtitle("Laplace") + geom_smooth(se=FALSE) +
  # scale_colour_ipsum("W") +
  ## scale_shape(name="Simulation V",solid=FALSE) +
  ylim(ymn2,ymx2) + ylab("log BF")+
  theme_cowplot(12) + 
  background_grid(major = "xy", minor = "none")+
  theme(plot.title = element_text(color="black", size=14, face="plain"))

# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  p3 + theme(legend.box="horizontal",
             legend.box.background = element_rect(colour = "black"),
             legend.box.margin = margin(6, 6, 6, 6))
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

ggsave("marg-lhood-bf.png",plot=final,height=8,width=10)

# Also plot the relationships
p.ls <- list()
k <- 1
for (i in ncol(cor.dat):1){
  for (j in 1:ncol(cor.dat)){
    v1 <- names(cor.dat)[i]
    v2 <- names(cor.dat)[j]
    p.ls[[k]] <- ggplot(cor.dat, aes_string(x=sym(v1), y=sym(v2), color=factor(cast.df$survey))) +
      geom_point() +
      geom_smooth(mapping=aes_string(x=sym(v1), y=sym(v2)), method='lm', inherit.aes = F) +
      ggtitle(paste0(v1,' vs ',v2)) +
      labs(color="Survey")
    k <- k + 1
  }
}
leg <- get_legend(p.ls[[1]])
for (k in 1:length(p.ls)){
  p.ls[[k]] <-  p.ls[[k]] + theme(legend.position = 'none')
}
for (k in c(5,9,10,13,14,15,17,18,19,20,21,22,23,24,25)){
  p.ls[[k]] <- geom_blank()
}
p.ls[[17]] <- leg
ggarrange(plotlist = p.ls, ncol=5, nrow=5)
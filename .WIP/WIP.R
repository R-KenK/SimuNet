# Tutorial draft  -----------------------------------------------------
Adj <- import_from_asnr("Mammalia","kangaroo_proximity_weighted",output = "adjacency",type = "upper")

total_scan <- 241
subset(SimuNet::asnr_network_df(),class == "Mammalia")
theo <- simu_scan(Adj,total_scan,scans.to.do = "all",mode = "upper")
theo
theo.sum <- summary(theo)
theo.sum
str(theo.sum)

plot_adj_cor(Adj[upper.tri(Adj)]/9000,theo.sum$theoretical.scaled[upper.tri(Adj)],
             main = "Association indices correlation\n(theoretical scans)",
             xlab = "Original association index",
             ylab = "Simulated association index"
)

plot_adj_cor <- function(X, Y,...) {
  plot(X,Y,...)
  abline(0,1,col = "red")
  text(0.7,0.1,labels = paste0("pearson's r = ",signif(cor.test(X,Y)$estimate,5)))
}

# unbiased random group-scan
para.halfgroup <- simu_samplingParam(Adj,total_scan,mode = "upper",scans.to.do = "all",group.scan_param = 0.5)
para.halfgroup
halfgroup <- simu_scan(sampling.param = para.halfgroup)
halfgroup
halfgroup.sum <- summary(halfgroup)
halfgroup.sum
str(halfgroup.sum)
hist(halfgroup.sum$group.sampled[upper.tri(halfgroup.sum$group.sampled)]/9000,
     breaks = 0.001*480:520,
     main = "edge sampling ratio",
     xlab = "proportion of time a dyad has been observed over the whole sampling effort"
)
abline(v = 0.5,col = "blue",lty = "dashed",lwd = 3)
text(0.51,17,labels = paste0("observation probability constant = ",max(para.halfgroup$obs.prob$P)))

plot_adj_cor(Adj[upper.tri(Adj)]/9000,halfgroup.sum$theoretical.scaled[upper.tri(Adj)],
             main = "Association indices correlation\n(theoretical scans)",
             xlab = "Original association index",
             ylab = "Simulated association index"
)

plot_adj_cor(halfgroup.sum$theoretical.scaled[upper.tri(Adj)],halfgroup.sum$group.scaled[upper.tri(Adj)],
             main = "Association indices correlation\n(group scans)",
             xlab = "Theoretical association index",
             ylab = "Empirical association index"
)

# unbiased random group-scan & even focal scan
para.both <- simu_samplingParam(Adj,total_scan,mode = "upper",scans.to.do = 1:3000,group.scan_param = .2,focal.scan_param = "even")
both <- simu_scan(sampling.param = para.both)
both
both.sum <- summary(both)
both.sum
plot_adj_cor(Adj[upper.tri(Adj)]/9000,
             both.sum$theoretical.scaled[upper.tri(Adj)])

plot_adj_cor(both.sum$theoretical.scaled[upper.tri(Adj)],
             both.sum$group.scaled[upper.tri(Adj)])
plot_adj_cor(both.sum$theoretical.scaled[upper.tri(Adj)],
             both.sum$focal.scaled[upper.tri(Adj)])


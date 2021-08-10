# Tutorial draft  -----------------------------------------------------
Adj <- import_from_asnr("Mammalia","bats_f",output = "adjacency",type = "upper")

total_scan <- 9000

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
para.both <- simu_samplingParam(Adj,total_scan,mode = "upper",
                                scans.to.do = 1:3000,group.scan_param = .2,focal.scan_param = "even")
both <- simu_scan(sampling.param = para.both)
both
both.sum <- summary(both)
both.sum
plot(both.sum)
str(both.sum)





#  triangular packed matrices vs regular vs vectors -----------------------

randomize_Adj <- function(n, total_scan,Adj.subfun = diag) {
  set.seed(42)
  nodes<- letters[1:n]
  Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
  Adj[non.diagonal(Adj)] <- round(runif(n*(n-1),0,total_scan))
  Adj[Adj.subfun(Adj)] <- NA
  Adj
}

bench.df <- rbind_lapply(
  c(10,30,50,100,150,200,250),
  function(n) {
    Adj <- randomize_Adj(n,42,function(Adj) {lower.tri(Adj,diag = TRUE)})
    Adj.packed <- pack_snPackMat(Adj,"upper")
    bench <- microbenchmark::microbenchmark(
      standard = {
        exp(Adj)
      },
      packed = {
        exp(Adj.packed)
      },
      unit = "ns",times = 50
    )
    df <- data.frame(n = n,
               method = bench$expr,
               type = "triangle",
               time = bench$time
    )
    Adj <- randomize_Adj(n,42,diag)
    Adj.packed <- pack_snPackMat(Adj,"directed")
    bench <- microbenchmark::microbenchmark(
      standard = {
        exp(Adj)
      },
      packed = {
        exp(Adj.packed)
      },
      unit = "ns",times = 50
    )
    rbind(df,
          data.frame(n = n,
                     method = bench$expr,
                     type = "non.diagonal",
                     time = bench$time
          )
    )
  }
)

library(ggplot2)

ggplot(bench.df,aes(n,time,colour = method))+
  facet_grid(type~.,scales = "free_y")+
  geom_smooth(aes(fill = method),method = "glm",formula = y ~ poly(x,degree = 2),alpha = .5)+
  geom_boxplot(aes(fill = method,group = interaction(n,method)),alpha = .7)+
  geom_point(fill = "white",shape = 21,alpha = .32)+
  scale_y_continuous(limits = c(0,NA))+
  scale_x_continuous(breaks = seq(25,from = 0, to = 300))+
  theme_bw()

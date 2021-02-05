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
plot_adj_cor(Adj[upper.tri(Adj)]/9000,
             both.sum$theoretical.scaled[upper.tri(Adj)])

plot_adj_cor(both.sum$theoretical.scaled[upper.tri(Adj)],
             both.sum$group.scaled[upper.tri(Adj)])
plot_adj_cor(both.sum$theoretical.scaled[upper.tri(Adj)],
             both.sum$focal.scaled[upper.tri(Adj)])


#  triangular packed matrices vs regular vs vectors -----------------------

Adj <- import_from_asnr("Mammalia","bats_f",output = "adjacency",type = "upper")
total_scan <- 9000


theo <- simu_scan(Adj,total_scan,scans.to.do = "all",mode = "upper")
theo
theo.sum <- summary(theo)
theo.sum
str(theo.sum)

para.both <- simu_samplingParam(Adj,total_scan,mode = "upper",
                                scans.to.do = "all",group.scan_param = .2,focal.scan_param = "even")

both <- simu_scan(sampling.param = para.both)
both
both.sum <- summary(both)
both.sum
str(both.sum)

bench.standard <- do.call(rbind,
                          lapply(
                            1:10,
                            function(i) {
                              cat(paste0("\r",i * 10,"%"))
                              start <- Sys.time()
                              theo <- simu_scan(Adj,total_scan,scans.to.do = "all",mode = "upper")
                              theo.sum <- summary(theo)
                              both <- simu_scan(sampling.param = para.both)
                              both.sum <- summary(both)
                              stop <- Sys.time()

                              data.frame(
                                i = i,
                                method = "standard",
                                # method = "packed",
                                time = difftime(stop,start,units = "secs"),
                                size.theo = format(object.size(theo),units = "Kb"),
                                size.theo.sum = format(object.size(theo.sum),units = "Kb"),
                                size.both = format(object.size(both),units = "Kb"),
                                size.both.sum = format(object.size(both.sum),units = "Kb")
                              )
                            }
                          )
)
#
# bench.packed <- do.call(rbind,
#                         lapply(
#                           1:10,
#                           function(i) {
#                             cat(paste0("\r",i * 10,"%"))
#                             start <- Sys.time()
#                             theo <- simu_scan(Adj,total_scan,scans.to.do = "all",mode = "upper")
#                             theo.sum <- summary(theo)
#                             both <- simu_scan(sampling.param = para.both)
#                             both.sum <- summary(both)
#                             stop <- Sys.time()
#
#                             data.frame(
#                               i = i,
#                               # method = "standard",
#                               method = "Matrix.packed",
#                               time = difftime(stop,start,units = "secs"),
#                               size.theo = format(object.size(theo),units = "Kb"),
#                               size.theo.sum = format(object.size(theo.sum),units = "Kb"),
#                               size.both = format(object.size(both),units = "Kb"),
#                               size.both.sum = format(object.size(both.sum),units = "Kb")
#                             )
#                           }
#                         )
# )

bench.homemade <- do.call(rbind,
                        lapply(
                          1:10,
                          function(i) {
                            cat(paste0("\r",i * 10,"%"))
                            start <- Sys.time()
                            theo <- simu_scan(Adj,total_scan,scans.to.do = "all",mode = "upper")
                            theo.sum <- summary(theo)
                            both <- simu_scan(sampling.param = para.both)
                            both.sum <- summary(both)
                            stop <- Sys.time()

                            data.frame(
                              i = i,
                              # method = "standard",
                              method = "SimuNet.packed",
                              time = difftime(stop,start,units = "secs"),
                              size.theo = format(object.size(theo),units = "Kb"),
                              size.theo.sum = format(object.size(theo.sum),units = "Kb"),
                              size.both = format(object.size(both),units = "Kb"),
                              size.both.sum = format(object.size(both.sum),units = "Kb")
                            )
                          }
                        )
)

bench <- rbind(bench.standard,bench.packed,bench.homemade)
bench$time <- as.numeric(bench$time)

bench$size.theo <- as.numeric(substr(bench$size.theo,1,nchar(bench$size.theo) - 3))
bench$size.theo.sum <- as.numeric(substr(bench$size.theo.sum,1,nchar(bench$size.theo.sum) - 3))
bench$size.both.sum <- as.numeric(substr(bench$size.both.sum,1,nchar(bench$size.both.sum) - 3))
bench$size.both <- as.numeric(substr(bench$size.both,1,nchar(bench$size.both) - 3))

library(ggplot2)

ggplot(bench,aes(method, as.numeric(time), fill = method))+
  geom_boxplot(colour = "black")+theme_bw()

ggplot(bench,aes(method, size.theo, fill = method))+
  geom_boxplot(colour = "black")+theme_bw()

median(bench[bench$method == "Matrix.packed",]$time) / median(bench[bench$method == "standard",]$time)
median(bench[bench$method == "Matrix.packed",]$size.theo) / median(bench[bench$method == "standard",]$size.theo)
median(bench[bench$method == "Matrix.packed",]$size.theo.sum) / median(bench[bench$method == "standard",]$size.theo.sum)
median(bench[bench$method == "Matrix.packed",]$size.both) / median(bench[bench$method == "standard",]$size.both)
median(bench[bench$method == "Matrix.packed",]$size.both.sum) / median(bench[bench$method == "standard",]$size.both.sum)

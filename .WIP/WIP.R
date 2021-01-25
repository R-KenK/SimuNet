# left over bits of code ----


minimize_NA<- function(scan){
  if(!is.empiScan(scan)) {scan}
  switch(scan$mode,
         "undirected" = , # as in `igraph`, consider this mode to be the same as `max`
         "max" = {
           resolvable_NA<- is.na(scan) & ((t(scan) == 1) & !is.na(t(scan)))  # scan[i,j] = NA & scan[j,i] = 1 =>  scan[i,j]<- 1
           ifelse(resolvable_NA,1,scan)
         },
         "min" = {
           resolvable_NA<- is.na(scan) & ((t(scan) == 0) & !is.na(t(scan)))  # scan[i,j] = NA & scan[j,i] = 0 =>  scan[i,j]<- 1
           ifelse(resolvable_NA,0,scan)
         },
         "plus" = ,
         "directed" = ,
         "upper" = ,
         "lower" =  ,
         "vector" = raw.scan
  )
}



# manual tests to see if `apply_mode` does what it should, including on weighted networks ---------------

# Not needed anymore because now `apply_mode` works on the `theoretical.scan.list` part of `scan` object to produce `empiScan` objects
# test.scan.NA<- matrix(c(0,1,NA,0,NA,
#                         NA,0,1,1,NA,
#                         0,1,0,0,1,
#                         0,0,1,0,0,
#                         NA,1,0,0,0),nrow = 5,byrow = TRUE,dimnames = list(letters[1:5],letters[1:5]))

# handcrafted example that contain all potential problematic symmetrization issue
test.raw.scan<- matrix(c(0,1,0,1,
                0,0,0,0,
                1,0,0,1,
                1,0,1,0),nrow = 4,byrow = TRUE,dimnames = list(letters[1:4],letters[1:4]))

apply_mode(test.raw.scan,"directed")
apply_mode(test.raw.scan,"undirected")
apply_mode(test.raw.scan,"max")
apply_mode(test.raw.scan,"min")
apply_mode(test.raw.scan,"plus")

# handcrafted example that contain all potential problematic symmetrization issue
test.raw.scan.weighted<- matrix(c(0,20,5,3,
                         10,0,2,0,
                         30,2,0,0,
                         3,0,1,0),nrow = 4,byrow = TRUE,dimnames = list(letters[1:4],letters[1:4]))

apply_mode(test.raw.scan.weighted,"directed")
apply_mode(test.raw.scan.weighted,"undirected")
apply_mode(test.raw.scan.weighted,"max")
apply_mode(test.raw.scan.weighted,"min")
apply_mode(test.raw.scan.weighted,"plus")


# Usage example (as of 0.3.0.9000) for single scans -----------------------

#' The user should mostly interact with `simu_X()` functions. As of now, the
#' user should first define sampling parameters of class `samplingParam` through
#' `simu_samplingParam()`, and then input them into `simu_scan()` to produce
#' simulated empirical scans.

set.seed(42)

n<- 5;nodes<- letters[1:n];total_scan<- 42;
Adj<- matrix(data = 0,nrow = n,ncol = n,dimnames = list(nodes,nodes))
Adj[non.diagonal(Adj)]<- sample(0:total_scan,n*(n-1),replace = TRUE)
Adj

# by default will simulate a directed theoretical scan
plot(simu_scan(Adj,total_scan))

# but other mode can be used:
simu_scan(Adj,total_scan,mode = "min")


# Users can generate sampling parameters through `simu_samplingParam` to use in
# `simu_scan`
para.group.constant<- simu_samplingParam(Adj,total_scan,mode =
                                         "max",group.scan_param = 0.42)
simu_scan(Adj,total_scan,sampling.param = para.group.constant)

# Users can also define functions to use trait- or network- based sampling
# biases for group-scan sampling (cf. ?simu_samplingParam)
obs.prob.trait.bias_fun<- function(i,j,Adj) {i+j} # comparable to a dyad-trait-based bias
para.group.trait.bias<- simu_samplingParam(Adj,total_scan,mode ="directed",
                                           group.scan_param = obs.prob.trait.bias_fun)
para.group.net.bias<- simu_samplingParam(Adj,total_scan,mode =
                                         "max",group.scan_param = function(i,j,Adj) {Adj*Adj})

simu_scan(Adj,total_scan,sampling.param = para.group.trait.bias)
simu_scan(Adj,total_scan,sampling.param = para.group.net.bias)

# or for biases regarding which focals to draw for focal-scan sampling (cf.
# ?simu_samplingParam)
focal.trait.bias_fun<- function(n,Adj) {1:n} # comparable to a dyad-trait-based bias
para.focal.trait.bias<- simu_samplingParam(Adj,
                                           total_scan,mode = "directed",
                                           focal.scan_param = focal.trait.bias_fun,
                                           scans.to.do = 3:20)
para.focal.net.bias<- simu_samplingParam(Adj,total_scan,mode = "max",
                                         focal.scan_param = function(n,Adj) {colSums(Adj*Adj)},
                                         scans.to.do = 20)

simu_scan(Adj,total_scan,sampling.param = para.focal.trait.bias)



simu_samplingParam(Adj,total_scan,focal.scan_param = "even")
simu_samplingParam(Adj,total_scan,mode = "max",group.scan_param = 0.42)
simu_samplingParam(Adj,total_scan,mode = "min",
                   group.scan_param = 0.42,
                   focal.scan_param = "random",scans.to.do = 1:4)


# Testing with Mamiko -----------------------------------------------------

#' TO WRITE
#'
#' @param path TO WRITE
#' @param output TO WRITE
#' @param type TO WRITE
#'
#' @return TO WRITE
#' @noRd
import_from_graphml<- function(path,output = c("graph","adjacency"),type = c("both", "upper", "lower")){
  output <- match.arg(output)
  G <- igraph::read_graph(path,format = "graphml")
  switch(output,
         "graph" = G,
         "adjacency" = {
           Adj <- igraph::get.adjacency(G,type = type,attr = "weight",sparse = FALSE);
           if(!is.null(igraph::vertex_attr(G,"name")) | !is.null(igraph::vertex_attr(G,"id"))) {
             if(is.null(igraph::vertex_attr(G,"name"))) {
               rownames(Adj) <- igraph::vertex_attr(G,"id")
             } else {
               rownames(Adj) <- igraph::vertex_attr(G,"name")
             }
             colnames(Adj) <- rownames(Adj)
           } else {
             rownames(Adj) <- as.character(1:nrow(Adj));colnames(Adj)<- as.character(1:ncol(Adj))
           }
           Adj
         }
  )
}

Adj <- import_from_graphml("C:/R/Git/asnr/Networks/Mammalia/bats_foodsharing_weighted/vampirebats_carter_mouth_licking_attribute.graphml",
                           "adjacency",
                           "upper")
source("C:/R/Git/asnr/Networks/Mammalia/bats_foodsharing_weighted/total_scan.R")
total_scan

theo <- simu_scan(Adj,total_scan,scans.to.do = "all",mode = "upper")
theo.sum <- summary(theo)
theo.sum$theoretical.scaled
plot(Adj,theo.sum$theoretical.sum)

# unbiased random group-scan
para.group <- simu_samplingParam(Adj,total_scan,mode = "upper",scans.to.do = "all",group.scan_param = "random")
group <- simu_scan(sampling.param = para.group)
group.sum <- summary(group)
plot(Adj,group.sum$theoretical.sum)
plot(group.sum$theoretical.scaled,group.sum$group.scaled)

# unbiased random group-scan & even focal scan
para.both <- simu_samplingParam(Adj,total_scan,mode = "upper",scans.to.do = "all",group.scan_param = "random",focal.scan_param = "even")
both <- simu_scan(sampling.param = para.both)
both.sum <- summary(both)
plot(Adj[upper.tri(Adj)],
     both.sum$theoretical.sum[upper.tri(both.sum$theoretical.sum)])

plot(both.sum$theoretical.scaled[upper.tri(both.sum$theoretical.scaled)],
     both.sum$group.scaled[upper.tri(both.sum$group.scaled)])
plot(both.sum$theoretical.scaled[upper.tri(both.sum$theoretical.scaled)],
     both.sum$focal.scaled[upper.tri(both.sum$focal.scaled)])

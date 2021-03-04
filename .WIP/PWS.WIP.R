# installation ------------------------------------------------------------
renv::install.packages("R-KenK/SimuNet@PWS")
devtools::install_github("R-KenK/SimuNet@PWS")
library(SimuNet)
library(ggplot2)

# importing data ------------------------------------------------------
# Own data
Adj.own <- matrix(c(0,90,60,56,1,
                40,0,0,20,3,
                20,0,0,10,2,
                80,0,3,0,5,
                0,0,0,0,0),
              ncol = 5,byrow = TRUE,
              dimnames = list(1:5,1:5)
)
total_scan.own <- 100
Adj.own

# importing from paper
obs.data <-
  c(  0,119, 40,13,38,16,18,37,
    119,  0,101, 0, 0,23, 0, 0,
     40,101,  0,31,74,94, 6,20,
     13,  0, 31, 0,39,78, 4,23,
     38,  0, 74,39, 0,33, 0,17,
     16, 23, 94,78,33, 0,10,33,
     18,  0,  6, 4, 0,10, 0,68,
     37,  0, 20,23,17,33,68, 0)
nodes_names <- c("A2","A","A1","B",
                 "B2","B1","C","C1")
Adj.massen <- matrix(obs.data,
              ncol = 8,
              byrow = TRUE,
              dimnames = list(nodes_names,
                              nodes_names
              )
)
total_scan.massen <- 1138
Adj.massen

# Importing from ASNR
## from URL
Adj.massen <- import_from_asnr(
  url = "https://github.com/bansallab/asnr/blob/master/Networks/Mammalia/rhesusmacaque_association_weighted/macaque_massen_contact_sits_attribute.graphml",
  output = "adj"
)
total_scan.massen <- 1138
Adj.massen

## from folder name
Adj.grant <- import_from_asnr("Mammalia",
                        "kangaroo_proximity_weighted",
                        output = "adjacency",type = "upper")
total_scan.grant <- 241 # number of observations reported by the author
# data retrieved from (Grant, 1973)
# DOI: 10.1016/S0003-3472(73)80004-1
Adj.grant

# Comparing theoretical and empirical scans -------------------------------

# Simulating theoretical scans
simulated.scans <- simu_scan(Adj.massen,total_scan.massen,scans.to.do = "all",mode = "upper")
simulated.scans
summary(simulated.scans)
plot(simulated.scans,method = "both",
     vertex.label = NA,vertex.size.min = 5,vertex.color = "#70c0b1", edge.color = "grey80")

# Simulating unbiased focal scans
## Setting sampling parameters: focal sampling
unbiased.focal <- simu_samplingParam(Adj.massen,total_scan.massen,"upper",focal.scan_param = "even",scans.to.do = "all")

## Running the empirical simulations: focal sampling
simulated.focals <- simu_scan(sampling.param = unbiased.focal)
simulated.focals
simulated.focals.summ <- summary(simulated.focals)
print(simulated.focals.summ,scaled = TRUE)

plot(simulated.focals,method = "both",
     vertex.label = NA,vertex.size.min = 6,vertex.color = "green")


## Setting sampling parameters: group & focal sampling
unbiased.both <- simu_samplingParam(Adj.grant,total_scan.grant,"upper",group.scan_param = 0.5,focal.scan_param = "even",scans.to.do = "all")

## Running the empirical simulations: group & focal sampling
simulated.both <- simu_scan(sampling.param = unbiased.both)
simulated.both
simulated.both.summ <- summary(simulated.both)
print(simulated.both.summ,scaled = TRUE)
plot(simulated.both,vertex.label = NA,vertex.size.min = 10,vertex.color = "lightblue")

# Get data ----------------------------------------------------------------
# theoretical vs focal

replicated.simu <- replicate(500,simu_scan(sampling.param = unbiased.both),simplify = FALSE)

globalCC.df.simple <- rbind_lapply(
  seq_along(replicated.simu),
  function(r) {
    simu <- replicated.simu[[r]]
    data.frame(
      replication = r,
      method = c("theoretical","group","focal"),
      GlobalCC = calculate_SNm(simu,
                               c("theoretical","group","focal"),
                               "scaled",
                               GlobalCC)
    )
  }
)

globalCC.df.simple$method <- factor(globalCC.df.simple$method,levels = c("theoretical","group","focal"))
ggplot(globalCC.df.simple,aes(method,GlobalCC,colour = method,fill = method))+
  geom_jitter(alpha = .3)+
  geom_boxplot(colour = "black",alpha = 0.4)+
  ylab("Global Clustering coefficient")+
  theme_bw()

strength.df.simple <- rbind_lapply(
  seq_along(replicated.simu),
  function(r) {
    simu <- replicated.simu[[r]]
    nodes <- rownames(simu$Adj)
    data.frame(
      replication = r,
      method = c("theoretical","group","focal"),
      node = factor(rep(nodes,3),levels = nodes),
      strength = calculate_SNm(simu,
                               c("theoretical","group","focal"),
                               "scaled",
                               compute.strength)
    )
  }
)

strength.df.simple$method <- factor(strength.df.simple$method,levels = c("theoretical","group","focal"))

ggplot(strength.df.simple,
       aes(node,strength,colour = method,fill = method,group = interaction(node,method)))+
  geom_jitter(alpha = .15)+
  geom_boxplot(colour = "black",alpha = 0.4)+
  xlab("individual")+ylab("node strength")+
  theme_bw()


# example of simulating bias ----------------------------------------------
central_bias.function <- function(i,j,Adj) {
  (compute.EV(Adj,"upper"))^2 # probabilty of observation proportional to the node's eigen-vector centrality power 5
}
biased.group <- simu_samplingParam(Adj.grant,total_scan.massen,"upper",
                                   group.scan_param = central_bias.function,
                                   scans.to.do = "all")

## Running the empirical simulations: biased group sampling
simulated.biased.group <- simu_scan(sampling.param = biased.group)
simulated.biased.group
simulated.biased.group.summ <- summary(simulated.biased.group)
print(simulated.biased.group.summ,scaled = TRUE)

plot(simulated.biased.group,method = "both",
     vertex.label = NA,vertex.size.min = 5,vertex.color = "tomato")

replicated.simu.biased <- replicate(500,simu_scan(sampling.param = biased.group),simplify = FALSE)
replicated.simu.unbiased <-
  replicate(500,
            {
              unbiased.group <- simu_samplingParam(Adj.grant,total_scan.massen,"upper",
                                                   group.scan_param = "random",
                                                   scans.to.do = "all")
              simu_scan(sampling.param = unbiased.group)
            },simplify = FALSE
  )

globalCC.df.bias <- rbind_lapply(
  1:500,
  function(r) {
    simu.biased <- replicated.simu.biased[[r]]
    simu.unbiased <- replicated.simu.unbiased[[r]]
    bias <- c("theoretical","biased","random")
    data.frame(
      replication = r,
      method = factor(bias,levels = bias),
      GlobalCC = c(
        calculate_SNm(simu.biased,
                      "theoretical",
                      "scaled",
                      GlobalCC),
        calculate_SNm(simu.biased,
                      "group",
                      "scaled",
                      GlobalCC),
        calculate_SNm(simu.unbiased,
                      "group",
                      "scaled",
                      GlobalCC)
      )
    )
  }
)

ggplot(globalCC.df.bias,aes(method,GlobalCC,colour = method,fill = method))+
  geom_jitter(alpha = .3)+
  geom_boxplot(colour = "black",alpha = 0.4)+
  scale_x_discrete(labels = c("all observed","central bias","random (unbiased)"))+
  ylab("Global Clustering coefficient")+xlab("Type of observation bias")+
  theme_bw()


# custom functions: node removal ------------------------------------------



# Simulating theoretical scans
replicated.node_removal <- replicate(500,
                                     simu_scan(Adj.grant,total_scan.massen,
                                               scans.to.do = "all",mode = "upper"),
                                     simplify = FALSE
)

EV.df <- rbind_lapply(
  seq_along(replicated.node_removal),
  function(r) {
    simu.theo <- summary(replicated.node_removal[[r]])
    simu.removed <- simu.theo
    simu.removed$theoretical.sum <-
      remove_central.node(simu.removed$theoretical.sum,
                          "upper",compute.EV
      )
    network <- c("theoretical","central removed")
    nodes <- c(rownames(simu.theo$theoretical.sum),
               rownames(simu.removed$theoretical.sum)
    )
    data.frame(
      replication = rep(r,33),
      network = factor(c(rep("theoretical",17),rep("central removed",16)),levels = network),
      node = factor(nodes,levels = rownames(simu.theo$theoretical.sum)),
      EV = c(
        {
          l <- compute.SNm.list(simu.theo,
                                "theoretical.sum",
                                "upper",
                                compute.EV)
          format_output("vector",l)
        },
        {
          l.bis <- compute.SNm.list(simu.removed,
                                "theoretical.sum",
                                "upper",
                                compute.EV)
          format_output("vector",l.bis)
        }
      )
    )
  }
)

ggplot(EV.df,
       aes(node,EV,colour = network,fill = network,group = interaction(node,network)))+
  geom_jitter(alpha = .15)+
  geom_boxplot(colour = "black",alpha = 0.4)+
  ylab("node eigen-vector centrality")+xlab("individual")+
  theme_bw()

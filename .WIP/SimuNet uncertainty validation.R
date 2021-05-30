# source library and custom functions -------------------------------------
rm(list = ls())
source(".WIP/validation_tools.R")

set.seed(42)

n <- 10L
N <- 20L
n.rep <- 350L
n.samp <- 1000L

# Generating P ---------------------------------------------------------------------------------
P <- generate_P.seq(n,mode = "directed")

## Extracting original P data ------------------------------------------------------------------
Pij.dt <-
  P %>%
  extract_xij(x.name = "p",mode = "directed",N = N,method = as.factor("original"))

P.node.dt <-
  P %>%
  get_original.node.metrics()

P.global.dt <-
  P %>%
  get_original.global.metrics()


# Generate inferred data for a single network observation --------------------------------------
A0.list <- generate_infered_networks(P,N,n.samp = n.samp,seed = 42)

## retrieve P_hat ----
P_hat.dt <-
  A0.list %>%
  extract_pijs.vec(method = c("GT","bbinom","SimuNet","boot"),n = n,N = N,n.samp = n.samp) %>%
  do.call(rbind,.)

P_hat.dt[,ij := factor(paste0(i,"-",j)),] %>%
  mutate(ij = factor(ij,levels = (Pij.dt[,ij:=paste0(i,"-",j),]$ij))) %>%
  ggplot(aes(ij,p,colour = method,fill = method))+
  geom_boxplot(aes(group = interaction(method,ij)),alpha = .5)+
  geom_point(data = Pij.dt,aes(ij,p),inherit.aes = FALSE,
             colour = "darkred",size = 3)+
  theme_minimal_vgrid()

## retrieve node metrics ----
node_metric.dt <-
  A0.list %>%
  extract_node.metrics.vec(method = c("GT","bbinom","SimuNet","boot"),
                           n = n,N = N,n.samp = n.samp) %>%
  do.call(rbind,.)

node_metric.dt %>%
  ggplot(aes(i,EV,colour = method,fill = method))+
  geom_boxplot(aes(group = interaction(method,i)),alpha = .5)+
  geom_point(data = P.node.dt,aes(i,EV),inherit.aes = FALSE,
             colour = "darkred",size = 3)+
  theme_cowplot()

node_metric.dt %>%
  ggplot(aes(i,CC,colour = method,fill = method))+
  geom_boxplot(aes(group = interaction(method,i)),alpha = .5)+
  geom_point(data = P.node.dt,aes(i,CC),inherit.aes = FALSE,
             colour = "darkred",size = 3)+
  theme_cowplot()

node_metric.dt %>%
  ggplot(aes(i,fbet,colour = method,fill = method))+
  geom_boxplot(aes(group = interaction(method,i)),alpha = .5)+
  geom_point(data = P.node.dt,aes(i,fbet),inherit.aes = FALSE,
             colour = "darkred",size = 3)+
  theme_cowplot()

## retrieve global metrics ----
global_metric.dt <-
  A0.list %>%
  extract_global.metrics.vec(method = c("GT","bbinom","SimuNet","boot"),
                             n = n,N = N,n.samp = n.samp) %>%
  do.call(rbind,.)

global_metric.dt %>%
  ggplot(aes(method,diam,colour = method,fill = method))+
  geom_boxplot(aes(group = method),alpha = .5)+
  geom_point(data = P.global.dt,aes(x = 2.5,diam),inherit.aes = FALSE,
             colour = "darkred",size = 3)+
  theme_cowplot()

global_metric.dt %>%
  ggplot(aes(method,GCC,colour = method,fill = method))+
  geom_boxplot(aes(group = method),alpha = .5)+
  geom_point(data = P.global.dt,aes(x = 2.5,diam),inherit.aes = FALSE,
             colour = "darkred",size = 3)+
  theme_cowplot()

## summarize data ----
node_metric.dt %>%
  melt_node.metrics() %>%
  calculate_node.CIs() %>%
  merge_with_original.nodes(P.node.dt) %>%
  compare_with_CI() %>%
  proportion_in_CI()

# Generate inferred data for several network observations --------------------------------------
## parallelization ----
library(snow)
cl <- snow::makeCluster(7)
snow::clusterEvalQ(cl,expr = {source(".WIP/validation_tools.R")})
snow::clusterExport(cl,list = list("n","N","n.rep","n.samp","P","Pij.dt","P.node.dt","P.global.dt"))

## Gather data ----
# repeated <-
#   infer_multiple_networks(P,N,n.rep,n.samp,cl = cl)

### across several N ----
N.list <- c(20,50,100,1000,2000)
repeated <-
  infer_across_N(P,N = N.list,n.rep,n.samp,cl = cl) %>%
  {
    list(
      P_hat.dt = rbind_lapply(.,function(r) r$P_hat.dt),
      node_metric.dt = rbind_lapply(.,function(r) r$node_metric.dt),
      global_metric.dt = rbind_lapply(.,function(r) r$global_metric.dt)
    )
  }

snow::stopCluster(cl)

saveRDS(repeated,".WIP/repeated.across.N.rds")
## proto data exploration ----
repeated %>%
  .$P_hat.dt %>%
  proportion_in_CI.P_hat() %>%
  {
    .[
      ,prop.in := in.CI / n.rep,
    ][
      ,prop.lower := lower.CI / n.rep,
    ][
      ,prop.higher := higher.CI / n.rep,
    ][]
  }

repeated %>%
  .$node_metric.dt %>%
  proportion_in_CI.nodes() %>%
  {
    .[
      ,prop.in := in.CI / n.rep,
    ][
      ,prop.lower := lower.CI / n.rep,
    ][
      ,prop.higher := higher.CI / n.rep,
    ][]
  }

repeated %>%
  .$global_metric.dt %>%
  proportion_in_CI.global() %>%
  {
    .[
      ,prop.in := in.CI / n.rep,
    ][
      ,prop.lower := lower.CI / n.rep,
    ][
      ,prop.higher := higher.CI / n.rep,
    ][]
  }

repeated %>%
  .$P_hat.dt %>%
  ggplot(aes(ij,p,colour = method,fill = method))+
  geom_boxplot(aes(group = interaction(method,ij)),alpha = .5)+
  geom_point(data = Pij.dt[,ij := factor(paste0(i,"-",j)),],aes(ij,p),inherit.aes = FALSE,
             colour = "darkred",size = 3)+
  theme_cowplot()

repeated %>%
  .$node_metric.dt %>%
  subset(metric == "str") %>%
  ggplot(aes(i,value,colour = method,fill = method))+
  geom_boxplot(aes(group = interaction(method,i)),alpha = .5)+
  geom_point(data = P.node.dt,aes(i,str),inherit.aes = FALSE,
             colour = "darkred",size = 3)+
  theme_cowplot()

repeated %>%
  .$node_metric.dt %>%
  subset(metric == "fbet") %>%
  ggplot(aes(i,value,colour = method,fill = method))+
  geom_boxplot(aes(group = interaction(method,i)),alpha = .5)+
  geom_point(data = P.node.dt,aes(i,fbet),inherit.aes = FALSE,
             colour = "darkred",size = 3)+
  theme_cowplot()

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
  infer_across_N(P,N = N.list,n.rep,n.samp,cl = cl,Pij.dt,P.node.dt,P.global.dt) %>%
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
# Generate inferred data for several network observations & n --------------------------------------
rm(list = ls())
source(".WIP/validation_tools.R")

set.seed(424242)
# n.list <- c(10,25,50)
n.list <- c(75,100)
N.list <- c(50,100,1000)
n.rep <- 105L
n.samp <- 100L

## parallelization ----
library(snow)
cl <- snow::makeCluster(7)
snow::clusterEvalQ(cl,expr = {source(".WIP/validation_tools.R")})

## Gather data ----
# repeated <-
#   infer_multiple_networks(P,N,n.rep,n.samp,cl = cl)

### across several N & n ----
repeated.n <-
  infer_across_n_N(n = n.list,N = N.list,n.rep,n.samp,cl = cl) %>%
  combine_inferred_across_n_N()

snow::stopCluster(cl)

saveRDS(repeated.n,".WIP/repeated.across.n.75.100.N.rds")

output.list <- c(".WIP/repeated.across.n.10.25.50.N.rds",".WIP/repeated.across.n.10.25.50.more.N.rds")
repeated.n <- readRDS_and_rbind(output.list)
repeated.n <- readRDS(".WIP/repeated.across.n.10.25.50.N.rds")

## proto data exploration ----
### calculate proportions ----
repeated.n %>%
  .$P_hat.dt %>%
  proportion_in_CI.P_hat()

repeated.n %>%
  .$node_metric.dt %>%
  proportion_in_CI.nodes()

repeated.n %>%
  .$global_metric.dt %>%
  proportion_in_CI.global()


### plotting ----
repeated.n %>%
  .$P_hat.dt %>%
  ggplot(aes(ij,p,colour = method,fill = method))+
  facet_grid(n~.)+
  geom_boxplot(aes(group = interaction(method,ij)),alpha = .5)+
  theme_cowplot()

repeated.n %>%
  .$node_metric.dt %>%
  subset(metric == "str") %>%
  ggplot(aes(i,value,colour = method,fill = method))+
  facet_grid(n~.)+
  geom_boxplot(aes(group = interaction(method,i)),alpha = .5)+
  theme_cowplot()

repeated.n %>%
  .$node_metric.dt %>%
  subset(metric == "fbet") %>%
  ggplot(aes(i,value,colour = method,fill = method))+
  facet_grid(n~.)+
  geom_boxplot(aes(group = interaction(method,i)),alpha = .5)+
  theme_cowplot()

repeated.n %>%
  .$P_hat.dt %>%
  proportion_in_CI.P_hat() %>%
  {
    .[
      ,by = .(n,N,method),
      .(
        p.in = median(prop.in),
        p.in.low = low_q(prop.in),
        p.in.up = up_q(prop.in),
        p.lower = median(prop.lower),
        p.lower.low = low_q(prop.lower),
        p.lower.up = up_q(prop.lower),
        p.higher = median(prop.higher),
        p.higher.low = low_q(prop.higher),
        p.higher.up = up_q(prop.higher)
      )
    ]
  } %>%
  ggplot(aes(N,p.in,colour = method,fill = method))+
  facet_grid(n~.)+
  geom_hline(yintercept = 0.95,lty = "dashed",colour = "darkred")+
  geom_ribbon(aes(ymin = p.in.low,ymax = p.in.up),colour = NA, alpha = 0.2,
              position = position_dodge(10))+
  geom_point(position = position_dodge(10))+
  scale_x_continuous(breaks = c(50,100,1000))+
  theme_minimal_grid()

repeated.n %>%
  .$P_hat.dt %>%
  proportion_in_CI.P_hat() %>%
  select(-c("i","j","in.CI","lower.CI","higher.CI","n.rep")) %>%
  melt.data.table(id.vars = c("ij","n","N","method"),measure.vars = c("prop.in","prop.lower","prop.higher"),variable.name = "prop") %>%
  {
    .[
      ,by = .(n,N,method,prop),
      .(
        value = mean(value),
        low = low_q(value),
        up = up_q(value)
      )
    ]
  } %>%
  ggplot(aes(N,value,colour = method,fill = method))+
  facet_grid(prop~n,scales = "free")+
  geom_hline(data = data.table(prop = c("prop.in","prop.lower","prop.higher"),
                               value = c(0.95,0.025,0.025)),
             aes(yintercept = value),lty = "dashed",colour = "darkred")+
  geom_ribbon(aes(ymin = low,ymax = up),colour = NA, alpha = 0.2,
              position = position_dodge(10))+
  geom_point(position = position_dodge(10))+
  scale_x_continuous(breaks = c(0,100,1000,250,500,750,1500,2000))+
  theme_half_open()

repeated.n %>%
  {.$P_hat.dt[,abs.err := abs(p - original),][]} %>%
  {
    .[,by = .(n,N,method),
      .(
        abs.err = mean(abs.err),
        low = low_q(abs.err),
        up = up_q(abs.err)
      )
    ]
  } %>%
  ggplot(aes(N,abs.err,colour = method,fill = method))+
  facet_grid(.~n,scales = "free")+
  geom_ribbon(aes(ymin = low,ymax = up),colour = NA, alpha = 0.1,
              position = position_dodge(10))+
  geom_point(position = position_dodge(10))+
  scale_x_continuous(breaks = c(0,100,1000,250,500,750,1500,2000))+
  theme_half_open()

repeated.n %>%
  {.$P_hat.dt[,abs.err := abs(p - original),][]} %>%
  subset(n == 25) %>%
  subset(ij %in% paste0(seq(1,25,len = 5),"-",seq(1,25,len = 5) + 1)) %>%
  {
    .[,by = .(n,N,method,ij),
      .(
        abs.err = mean(abs.err),
        low = low_q(abs.err),
        up = up_q(abs.err)
      )
    ]
  } %>%
  ggplot(aes(N,abs.err,colour = method,fill = method))+
  facet_grid(.~ij,scales = "free")+
  geom_ribbon(aes(ymin = low,ymax = up),colour = NA, alpha = 0.1,
              position = position_dodge(10))+
  geom_point(position = position_dodge(10))+
  scale_x_continuous(breaks = c(0,100,1000,250,500,750,1500,2000))+
  theme_half_open()


repeated.n %>%
  .$node_metric.dt %>%
  proportion_in_CI.nodes() %>%
  select(-c("in.CI","lower.CI","higher.CI","n.rep")) %>%
  melt.data.table(id.vars = c("i","n","N","method","metric"),
                  measure.vars = c("prop.in","prop.lower","prop.higher"),
                  variable.name = "prop") %>%
  {
    .[
      ,by = .(n,N,method,metric,prop),
      .(
        value = mean(value),
        low = low_q(value),
        up = up_q(value)
      )
    ]
  } %>%
  subset(metric == "bet") %>%
  ggplot(aes(N,value,colour = method,fill = method))+
  facet_grid(prop~n,scales = "free")+
  geom_hline(data = data.table(prop = c("prop.in","prop.lower","prop.higher"),
                               value = c(0.95,0.025,0.025)),
             aes(yintercept = value),lty = "dashed",colour = "darkred")+
  geom_ribbon(aes(ymin = low,ymax = up),colour = NA, alpha = 0.2,
              position = position_dodge(10))+
  geom_point(position = position_dodge(10))+
  scale_x_continuous(breaks = c(0,100,1000,250,500,750,1500,2000))+
  theme_half_open()


repeated.n %>%
  .$node_metric.dt %>%
  subset(metric == "str" & n == 25 & i %in% seq(1,25,len = 5)) %>%
  {.[,abs.err := abs(value - original),][]} %>%
  {
    .[,by = .(n,N,method,i),
      .(
        abs.err = mean(abs.err),
        low = low_q(abs.err),
        up = up_q(abs.err)
      )
    ]
  } %>%
  ggplot(aes(N,abs.err,colour = method,fill = method))+
  facet_grid(.~i,scales = "free")+
  geom_ribbon(aes(ymin = low,ymax = up),alpha = 0.1,colour = NA,
              position = position_dodge(0))+
  geom_point(position = position_dodge(0))+
  scale_x_continuous(breaks = c(0,100,1000,250,500,750,1500,2000))+
  theme_half_open()

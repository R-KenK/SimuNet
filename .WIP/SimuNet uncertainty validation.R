# source library and custom functions -------------------------------------
rm(list = ls())
source(".WIP/validation_tools.R")

set.seed(42)

n <- 4L
N <- 20L
n.rep <- 10L
n.samp <- 100L

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

# parallelization
cl <- makeCluster(7)
snow::clusterEvalQ(cl,expr = {source(".WIP/validation_tools.R")})
snow::clusterExport(cl,list = list("n","N","n.rep","n.samp","P","Pij.dt","P.node.dt","P.global.dt"))

stopCluster(cl)


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
repeated <-
  infer_multiple_networks(P,N,n.rep,n.samp,cl = cl)

repeated %>%
  rbind_lapply(function(r) r$P_hat.dt) %>%
  proportion_in_CI.P_hat()

repeated %>%
  rbind_lapply(function(r) r$node_metric.dt) %>%
  proportion_in_CI.nodes()

repeated %>%
  rbind_lapply(function(r) r$global_metric.dt) %>%
  proportion_in_CI.global()


repeated %>%
  rbind_lapply(function(r) r$node_metric.dt) %>%
  subset(metric == "str") %>%
  ggplot(aes(i,value,colour = method,fill = method))+
  geom_boxplot(aes(group = interaction(method,i)),alpha = .5)+
  geom_point(data = P.node.dt,aes(i,str),inherit.aes = FALSE,
             colour = "darkred",size = 3)+
  theme_cowplot()

repeated %>%
  rbind_lapply(function(r) r$P_hat.dt) %>%
  proportion_in_CI.P_hat()

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
  proportion_in_CI(by = .(i,method,metric))





  {
    .[,by = .(n,N,method,i),
      .(EV = median(EV),EV.low = low_q(EV),EV.up = up_q(EV),
        CC = median(CC),CC.low = low_q(CC),CC.up = up_q(CC),
        str = median(str),str.low = low_q(str),str.up = up_q(str),
        bet = median(bet),bet.low = low_q(bet),bet.up = up_q(bet),
        fbet = median(fbet),fbet.low = low_q(fbet),fbet.up = up_q(fbet),
        n.samp = .N
      )
    ]
  } %>%
  ggplot(aes(method,CC,colour = method,fill = method))+
  facet_grid(.~i)+
  geom_linerange(aes(ymin = CC.low,ymax = CC.up))+
  geom_point()+
  theme_minimal_grid()









# old draft ----
library(snow)
library(pbapply)

cl <- makeCluster(7)
snow::clusterExport(cl,list = ls())
snow::clusterEvalQ(cl,expr = {library(data.table);library(dplyr);library(ggplot2);library(cowplot);library(extraDistr)})


CI.test <-
  pbreplicate(n = 105,simplify = FALSE,run_experiment(P,N,30,100),cl = cl) %>% do.call(rbind,.) %>%
  cbind(r2 = rep(1:14,each = (n^2 - n) * 2 * 30 * 100))
stopCluster(cl)

CI.test <- readRDS(".WIP/CI.test.rds")
low.bbinom <- function(A,N,CI.level = 0.95,alpha.prior = 0.5,beta.prior = 0.5,N.new = N) {
  extraDistr::pbbinom(0:N.new,N.new,alpha = A + alpha.prior,beta = N - A + beta.prior) %>% {
    max(which(. <= (1 - CI.level) / 2) - 1L,0L)
  }
}


up.bbinom <- function(A,N,CI.level = 0.95,alpha.prior = 0.5,beta.prior = 0.5,N.new = N) {
  if (A == N) {return(N.new)}
  extraDistr::pbbinom(0:N.new,N.new,alpha = A + alpha.prior,beta = N - A + beta.prior) %>%
    {max(which(. <= 1 - (1 - CI.level) / 2) - 1L)}
}

ci.bbinom <- function(A,N,CI.level = 0.95,alpha.prior = 0.5,beta.prior = 0.5,N.new = N) {
  extraDistr::pbbinom(0:N.new,N.new,alpha = A + alpha.prior,beta = N - A + beta.prior) %>%
    {c(max(which(. <= (1 - CI.level) / 2) - 1L,0L),
       ifelse(A == N,N.new,max(which(. <= 1 - (1 - CI.level) / 2) - 1L))
       )}
}

low.bbinomv <- Vectorize(low.bbinom,vectorize.args = "A")
up.bbinomv <- Vectorize(up.bbinom,vectorize.args = "A")
ci.bbinomv <- Vectorize(ci.bbinom,vectorize.args = "A")

CI.test[r2 == 1 & method == "Beta-binomial"
        ][,c("p.lower.beta","p.upper.beta") := .(low.bbinomv(a,N),up.bbinomv(a,N)),]


CI.test2 <- CI.test
CI.test2[,c("p.lower.beta","p.upper.beta") := .(qbeta(0.025,a +.5,N - a + .5),qbeta(0.975,a +.5,N - a + .5)),]

CI.test3 <- CI.test2[r2 == 1 & method == "Beta-binomial"][,in.CI := ifelse(p.lower.beta <= p & p <= p.upper.beta,1,0),]
CI.test3[,.(in.CI = sum(in.CI)/.N,n = .N),by = .(i,j,p,method,N,r1)] %>% {
  .[,.(in.CI = median(in.CI),lower = quantile(in.CI,.025),upper = quantile(in.CI,.975)),by = .(i,j,p,method,N)]
}

CI.stats <- CI.test[,by = .(r2,i,j,p,method,N,r1),.(p.lower = low(a / N),p.upper = up(a / N))]


CI.stats[,in.CI := ifelse(p.lower <= p & p <= p.upper,1,0),][]
CI.in <- CI.stats[,by = .(r2,i,j,p,method,N),.(prop.in.CI = sum(in.CI) / .N)]
CI.in[,by = .(i,j,p,method,N),.(prop.in.CI = median(prop.in.CI),lower = quantile(prop.in.CI,.025),upper = quantile(prop.in.CI,.975))]

test.dt[,by = .(i,j,p,method,N,r),.(p.lower = quantile(a / N,0.025),p.upper = quantile(a / N,0.975))] %>%
  {.[,in.CI := ifelse(p.lower <= p & p <= p.upper,1,0),][]} %>%
  {.[,by = .(i,j,p,method,N),.(prop.in.CI = sum(in.CI) / .N)]}

test.dt[,by = .(i,j,p,method,N,r),.(p.lower = quantile(a / N,0.025),p.upper = quantile(a / N,0.975))] %>%
  {.[,lower.CI := ifelse(p < p.lower,1,0),][]} %>%
  {.[,by = .(i,j,p,method,N),.(prop.lower.CI = sum(lower.CI) / .N)]}

subset(i == 4 & j == 3 & method == "Beta-binomial" & r == 1)

test.dt %>%
  subset(i == 4 & j == 3 & method == "Beta-binomial" & r == 1) %>%
  {.$a / .$N} %>%
  quantile(.,c(0.025,0.975))

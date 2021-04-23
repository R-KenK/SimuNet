# loading required tools and packages -----------------------------------------------
source(".WIP/P.obs.validation_tools.R")

# Simulation validation ---------------------------------------------------


## Generating "black-boxed" variables -------------------------------------

### toy-example variables -------------------------------------------------
set.seed(42)
n <- 5L
N <- 100L
n.rep <- 500L
node_names <- as.character(1:n)

### Drawing probability matrix of association -----------------------------
P0 <- draw_P0(n = n,method = "seq",min = 0,max = 1,node_names = node_names)
P0


### Generating Ground-truth adjacency matrix distribution -----------------
A.GT.list <-
  replicate(n = n.rep,
            simplify = FALSE,
            draw_Adj_from_binom(N,P0)
  )

## Generating researcher's empirical data ---------------------------------

### Drawing an empirical binary scan list ---------------------------------
scanList0 <- draw_scanList(P0,N = N)
scanList0

### Sum them into weighted adjacency matrix -------------------------------
A0 <- sum_scanList(scanList0)
A0


## Simulate scan lists and adjacency matrices from observed data ----------


### Simulations from A0 ---------------------------------------------------


#### Hierarchical model (beta probability, binomial for each scan) --------
A.sim.scanList <-
  replicate(n = n.rep,
            simplify = FALSE,
            A0 %>% draw_scanList(A.obs = .,N = N) %>% sum_scanList()
  )


#### Beta binomial model --------------------------------------------------
A.sim.bbinom <-
  replicate(n = n.rep,
            simplify = FALSE,
            A0 %>% draw_Adj_from_betabinom(A.obs = ,N = N)
  )


### Bootstrapping from scanList0 ------------------------------------------
A.bootstrap <-
  replicate(
    n = n.rep,
    simplify = FALSE,
    scanList0 %>% resample_scanList() %>% sum_scanList()
  )

## Compare GT and simulated data ------------------------------------------

### Extract data ----------------------------------------------------------
P0.df <- extract_xij(P0,x.name = "p")
A0.df <- extract_xij(A0,N = N,x.name = "a0",n.a0 = 1)

Comp.aij <-
  rbind(
    format_and_merge_xij.df(A.GT.list,N,"GT",A0.df),
    format_and_merge_xij.df(A.sim.scanList,N,"SimuNet",A0.df),
    format_and_merge_xij.df(A.sim.bbinom,N,"beta binomial",A0.df),
    format_and_merge_xij.df(A.bootstrap,N,"bootstrap",A0.df)
  )

Comp.aij$method <-
  Comp.aij$method %>%
  factor(levels = c("GT","beta binomial","SimuNet","bootstrap"))


### Plot data -------------------------------------------------------------


bw <- 1 / (3 * n * (n - 1))

A0.df <-
  Comp.aij %>% data.table %>%  {
    .[,by = .(N,method,i,j),
      .(a0 = mean(a0),
        type = as.factor(ifelse(method == "GT","GT","empirical")))
    ]
  }

CIs <-
  Comp.aij %>%
  group_by(N,n.a0,method,i,j,a0) %>%
  mutate(lower = quantile(a,0.025),
            upper = quantile(a,0.975))

Comp.aij %>%
  ggplot(aes(a / N,method,colour = method,fill = method))+
  facet_grid(i~j)+
  geom_density_ridges(stat = "binline",
                      binwidth = bw,
                      alpha = 0.4,
                      scale = 1.8)+
  geom_hline(yintercept = 1,colour = "black")+
  scale_x_continuous(limits = c(0 - bw,1 + bw),
                     breaks = c(0,0.5,1),labels = c("0","0.5","1"))+
  scale_y_discrete(expand = expansion(mult = c(0,0.5)))+
  xlab("a / N")+
  geom_vline(data = P0.df,aes(xintercept = p),
             lty = "dashed",colour = "tomato")+
  geom_point(data = A0.df,aes(x = a0 / N),
             colour = "royalblue")+
  geom_errorbarh(data = CIs,aes(xmin = lower / N,xmax = upper / N),height = 0.2,colour = "royalblue")+
  cowplot::theme_half_open()+
  theme(panel.grid.major.x = element_line(colour = "grey90"))

# Simulating across N -----------------------------------------------------
rm(list = ls())
source(".WIP/P.obs.validation_tools.R")

set.seed(42)
n <- 4L
n.rep <- 100L
node_names <- as.character(1:n)

N.to.do <- c(10,20,50,100,200,500)
# N.to.do <- c(50)
n.a0.to.do <- 1:1000
parameters <- expand.grid(N = N.to.do,n.a0 = n.a0.to.do)

P0 <- draw_P0(n = n,method = "seq",min = 0,max = 1,node_names = node_names)
P0.df <- extract_xij(P0,x.name = "p")

library(pbapply)
library(snow)

cl <- snow::makeCluster(7)
snow::clusterExport(cl,ls())

clusterEvalQ(cl,expr = {library(data.table);library(dplyr);library(extraDistr)})

Comp.aij.N <-
  pblapply(
    1:nrow(parameters),
    function(para) {
      N <- parameters$N[para];n.a0 <- parameters$n.a0[para]
      cat("\rN =",N,"- n.a0 =",n.a0,"             ")
      scanList0 <- draw_scanList(P0,N = N)
      A.GT.list <-
        replicate(n = n.rep,
                  simplify = FALSE,
                  draw_Adj_from_binom(N,P0)
        )
      A0 <- sum_scanList(scanList0)
      A.sim.scanList <-
        replicate(n = n.rep,
                  simplify = FALSE,
                  A0 %>% draw_scanList(A.obs = .,N = N) %>% sum_scanList()
        )
      A.sim.bbinom <-
        replicate(n = n.rep,
                  simplify = FALSE,
                  A0 %>% draw_Adj_from_betabinom(A.obs = ,N = N)
        )
      A.bootstrap <-
        replicate(
          n = n.rep,
          simplify = FALSE,
          scanList0 %>% resample_scanList() %>% sum_scanList()
        )
      A0.df <- extract_xij(A0,N = N,x.name = "a0",n.a0 = n.a0)

      Comp.aij <-
        rbind(
          format_and_merge_xij.df(A.GT.list,N,"GT",A0.df),
          format_and_merge_xij.df(A.sim.scanList,N,"SimuNet",A0.df),
          format_and_merge_xij.df(A.sim.bbinom,N,"beta binomial",A0.df),
          format_and_merge_xij.df(A.bootstrap,N,"bootstrap",A0.df)
        )

      Comp.aij$method <-
        Comp.aij$method %>%
        factor(levels = c("GT","beta binomial","SimuNet","bootstrap"))
      data.table(Comp.aij)
    },cl = cl
  )

snow::stopCluster(cl)
Comp.aij.N %>% setDT
Comp.aij.N <-
  Comp.aij.N %>%
  do.call(rbind,.)

a0.N.df <-
  Comp.aij.N %>% setDT %>%  {
    .[n.a0 == 1,by = .(N,method,i,j),
      .(a0 = mean(a0),
      type = as.factor(ifelse(method == "GT","GT","empirical")))
      ]
  }


bw <- 1 / (3 * n * (n - 1))

Comp.aij.N %>%
  subset(N %in% c(20,50,500) & n.a0 == 1) %>%
  mutate(type = as.factor(ifelse(method == "GT","GT","empirical"))) %>%
  mutate(type = relevel(type,"GT")) %>%
  group_by(N,n.a0,method,i,j,a0) %>%
  mutate(lower = quantile(a,0.025),
         upper = quantile(a,0.975)) %>%

  ggplot(aes(a / N,interaction(type,N),colour = method,fill = method))+
  facet_grid(i~j,scale = "free")+
  geom_density_ridges(stat = "binline",
                      binwidth = 1 / 50,
                      position = position_dodge(2 * bw),
                      alpha = 0.1,
                      scale = 1.2)+
  # geom_histogram(aes(y = stat(count) / sum(count)),
  #                alpha = 0.2,
  #                colour = NA,
  #                position = position_dodge(0),
  #                binwidth = 1 / 20)+
  geom_hline(yintercept = 1,colour = "black")+
  scale_x_continuous(limits = c(0 - bw,1 + bw),
                     breaks = c(0,0.5,1),labels = c("0","0.5","1"))+
  scale_y_discrete(expand = expansion(mult = c(0,0.5)),
                   labels = rep(c("N = 20","N = 50","N = 500"),each = 2))+
  geom_vline(data = P0.df,aes(xintercept = p),
             lty = "dashed",colour = "tomato")+
  geom_point(data = a0.N.df[N %in% c(20,50,500)],aes(x = a0 / N),
             colour = "royalblue")+
  geom_errorbarh(aes(xmin = lower / N,xmax = upper / N),height = 0.2)+
  # geom_label(data = P0.df,aes(x = 0.5,y = 7,label = paste0("p = ",round(p,2))))+
  xlab("a / N")+
  theme_ridges(center_axis_labels = TRUE,grid = FALSE)+
  # cowplot::theme_half_open()+
  theme(panel.grid.major.x = element_line(colour = "grey90"))

Comp.aij.N %>%
  subset(N %in% c(20,50,500)) %>%
  merge(P0.df,by = c("i","j")) %>%
  group_by(N,n.a0,method,i,j) %>%
  mutate(lower = quantile(a,0.025),
         upper = quantile(a,0.975),
         in.CI = ifelse(p >= lower / N & p <= upper / N,1,0)) %>%
  group_by(N,method,i,j,p) %>%
  summarize(p.value = sum(in.CI) / n()) %>%
  # group_by(N,method) %>%
  # summarise(p.value = mean(p.value),
  #           lower = quantile(p.value,0.025),
  #           upper = quantile(p.value,0.975)
  # ) %>%
  data.table %>%
  ggplot(aes(method,p.value))



a0.N.summ <-
  Comp.aij.N %>% setDT %>%  {
    .[,by = .(N,method,i,j),
      .(a0 = mean(a0),
        type = as.factor(ifelse(method == "GT","GT","empirical")))
    ]
  }

Comp.aij.N %>%
  subset(N %in% c(20,50,500)) %>%
  mutate(type = as.factor(ifelse(method == "GT","GT","empirical"))) %>%
  mutate(type = relevel(type,"GT")) %>%
  ggplot(aes(a / N,interaction(type,N),colour = method,fill = method))+
  facet_grid(i~j,scale = "free")+
  geom_density_ridges(stat = "binline",
                      binwidth = 1 / 50,
                      position = position_dodge(2 * bw),
                      alpha = 0.1,
                      scale = 1.2)+
  # geom_histogram(aes(y = stat(count) / sum(count)),
  #                alpha = 0.2,
  #                colour = NA,
  #                position = position_dodge(0),
  #                binwidth = 1 / 20)+
  geom_hline(yintercept = 1,colour = "black")+
  scale_x_continuous(limits = c(0 - bw,1 + bw),
                     breaks = c(0,0.5,1),labels = c("0","0.5","1"))+
  scale_y_discrete(expand = expansion(mult = c(0,0.5)),
                   labels = rep(c("N = 20","N = 50","N = 500"),each = 2))+
  geom_vline(data = P0.df,aes(xintercept = p),
             lty = "dashed",colour = "tomato")+
  geom_point(data = a0.N.summ[N %in% c(20,50,500)],aes(x = a0 / N),
             colour = "royalblue")+
  # geom_label(data = P0.df,aes(x = 0.5,y = 7,label = paste0("p = ",round(p,2))))+
  xlab("a / N")+
  cowplot::theme_half_open()+
  theme(panel.grid.major.x = element_line(colour = "grey90"))



# Old draft ---------------------------------------------------------------
#
#
# set.seed(42)
# n <- 5L
# node_names <- as.character(1:n)
# # P.ori <- draw_P.ori(n,rng.vec = rnorm(n,50,10),node_names = node_names,min = 0,max = .001,noise.sd = 0.2,adjust.min = 0)
# P.ori <- matrix(0,ncol = n,nrow = n,dimnames = list(node_names,node_names)) %>% {diag(.) <- 0; non.diagonal(.) <- seq(0,1,length.out = n * (n - 1));.}
#
# # compare methods with overall stochasticity (new A.obs at each rep)
# test.df <- rbind_lapply(
#   c(10,25,50,75,100),
#   function(N) {
#     node_names <- as.character(1:n)
#
#     cat("\nN = ",N,"\n")
#     rbind_lapply(
#       1:250,
#       function(r) {
#         A.obs <- rbinom_p(N,P.ori)
#
#         A.sim.simple <-
#           simulate_scan_old(A.obs,N) %>%
#           {Reduce("+",get_scanList(.,"simple"))}
#
#         A.sim.uniform <-
#           simulate_scan_old(A.obs,N) %>%
#           {Reduce("+",get_scanList(.,"uniform"))}
#
#         A.sim.SimuNet <-
#           simulate_scan_old(A.obs,N) %>%
#           {Reduce("+",get_scanList(.,"SimuNet"))}
#
#         A.sim.new <-
#           simulate_scan(A.obs,N) %>%
#           {Reduce("+",get_scanList(.,"A"))}
#
#         cat("\r",r)
#         expand.grid(i = 1:n,j = 1:n) %>% {
#           .$N <- N
#           .$theo <- P.ori[cbind(.$i,.$j)]
#           .$observed <- A.obs[cbind(.$i,.$j)]
#           .$uniform <- A.sim.uniform[cbind(.$i,.$j)]
#           .$simple <- A.sim.simple[cbind(.$i,.$j)]
#           .$SimuNet <- A.sim.SimuNet[cbind(.$i,.$j)]
#           .$bayesian <- A.sim.new[cbind(.$i,.$j)]
#           .
#         } %>%
#           data.table(
#             rep = r,
#             .
#           )
#       }
#     )
#   }
# )
#
# test.df %>%
#   subset(i != j) %>%
#   melt(id.vars = 1:6,measure.vars = 7:10,variable.name = "method",value.name = "simulated") %>%
#   .[,.(simulated = mean(simulated),
#        lower.sim = quantile(simulated,probs = 0.025),
#        upper.sim = quantile(simulated,probs = 0.975),
#        observed = mean(observed),
#        lower = quantile(observed,probs = 0.025),
#        upper = quantile(observed,probs = 0.975)),by = .(i,j,N,theo,method)] %>%
#   # ggplot(aes(observed / N,simulated / N,colour = as.factor(N)))+
#   ggplot(aes(theo,simulated / N,colour = as.factor(N)))+
#   # ggplot(aes(theo,observed / N,colour = as.factor(N)))+
#   facet_grid(method~N)+
#   geom_hline(yintercept = c(0,1))+
#   geom_vline(xintercept = 0)+
#   geom_abline(intercept = 0,slope = 1,colour = "grey50",lty = "dashed")+
#   geom_hline(aes(yintercept = .5 / (N + 1)),lty = "dashed",colour = "blue")+
#   geom_hline(aes(yintercept = (N + .5) / (N + 1)),lty = "dashed",colour = "blue")+
#   geom_hline(aes(yintercept = 1 / N),lty = "dotted",colour = "darkred")+
#   geom_hline(aes(yintercept = 1 - 1 / N),lty = "dotted",colour = "darkred")+
#   geom_smooth(formula = y ~ x,aes(fill = as.factor(N)),method = "lm",alpha = 0.3)+
#   geom_linerange(aes(ymin = lower.sim / N,ymax = upper.sim / N),alpha = 0.5,colour = "grey70")+
#   geom_point(alpha = 0.4)+
#   scale_y_continuous(limits = c(0,1),expand = expansion(mult = c(0,0.01)))+
#   scale_x_continuous(limits = c(0,1),expand = expansion(mult = c(0,0.01)))+
#   cowplot::theme_half_open()
#
# # compare methods with simulation stochasticity (a single A.obs drawn)
# set.seed(42)
# n <- 5L
# node_names <- as.character(1:n)
# # P.ori <- draw_P.ori(n,rng.vec = rnorm(n,50,10),node_names = node_names,min = 0,max = .001,noise.sd = 0.2,adjust.min = 0)
# P.ori <-
#
# test.df <- rbind_lapply(
#   c(10,25,50,75,100),
#   function(N) {
#     node_names <- as.character(1:n)
#
#     A.obs <- rbinom_p(N,P.ori)
#
#     cat("\nN = ",N,"\n")
#     rbind_lapply(
#       1:250,
#       function(r) {
#
#         A.sim.old <-
#           simulate_scan_old(A.obs,N) %>% {
#             list(
#               uniform = Reduce("+",get_scanList(.,"uniform")),
#               simple = Reduce("+",get_scanList(.,"simple")),
#               SimuNet = Reduce("+",get_scanList(.,"SimuNet"))
#             )
#           }
#
#         A.sim.uniform <- A.sim.old$uniform
#         A.sim.simple <- A.sim.old$simple
#         A.sim.SimuNet <- A.sim.old$SimuNet
#
#         A.sim.new <-
#           simulate_scan(A.obs,N) %>%
#           {Reduce("+",get_scanList(.,"A"))}
#
#         cat("\r",r)
#         expand.grid(i = 1:n,j = 1:n) %>% {
#           .$N <- N
#           .$theo <- P.ori[cbind(.$i,.$j)]
#           .$observed <- A.obs[cbind(.$i,.$j)]
#           .$uniform <- A.sim.uniform[cbind(.$i,.$j)]
#           .$simple <- A.sim.simple[cbind(.$i,.$j)]
#           .$SimuNet <- A.sim.SimuNet[cbind(.$i,.$j)]
#           .$bayesian <- A.sim.new[cbind(.$i,.$j)]
#           .
#         } %>%
#           data.table(
#             rep = r,
#             .
#           )
#       }
#     )
#   }
# )
#
# test.df %>%
#   subset(i != j) %>%
#   melt(id.vars = 1:6,measure.vars = 7:10,variable.name = "method",value.name = "simulated") %>%
#   .[,by = .(i,j,N,theo,method),
#     .(simulated = mean(simulated),
#       lower.sim = quantile(simulated,probs = 0.025),
#       upper.sim = quantile(simulated,probs = 0.975),
#       observed = mean(observed),
#       lower = quantile(observed,probs = 0.025),
#       upper = quantile(observed,probs = 0.975))] %>%
#   ggplot(aes(observed / N,simulated / N,colour = as.factor(N)))+
#   # ggplot(aes(theo,simulated / N,colour = as.factor(N)))+
#   # ggplot(aes(theo,observed / N,colour = as.factor(N)))+
#   facet_grid(method~N)+
#   geom_hline(yintercept = c(0,1))+
#   geom_vline(xintercept = 0)+
#   geom_abline(intercept = 0,slope = 1,colour = "grey50",lty = "dashed")+
#   geom_hline(aes(yintercept = .5 / (N + 1)),lty = "dashed",colour = "blue")+
#   geom_hline(aes(yintercept = (N + .5) / (N + 1)),lty = "dashed",colour = "blue")+
#   geom_hline(aes(yintercept = 1 / N),lty = "dotted",colour = "darkred")+
#   geom_hline(aes(yintercept = 1 - 1 / N),lty = "dotted",colour = "darkred")+
#   geom_smooth(formula = y ~ x,aes(fill = as.factor(N)),method = "lm",alpha = 0.3)+
#   geom_linerange(aes(ymin = lower.sim / N,ymax = upper.sim / N),alpha = 0.5,colour = "grey70")+
#   geom_point(alpha = 0.4)+
#   scale_y_continuous(limits = c(0,1),expand = expansion(mult = c(0,0.01)))+
#   scale_x_continuous(limits = c(0,1),expand = expansion(mult = c(0,0.01)))+
#   cowplot::theme_half_open()
#
#
#
# N <- 10L
# X <- 0:N
#
# diff.df <-
#   rbind_lapply(
#   c("uniform","simple","SimuNet","bayesian"),
#   function(method) {
#     cat("\nmethod:",method)
#     rbind_lapply(
#       seq(0,1,by = 0.1),
#       function(p) {
#         cat("\np = ",p,"\n")
#         rbind_lapply(
#           1:100,
#           function(r) {
#             cat("\r",r)
#             draw_p.hat(X,N,method) %>%
#               calculate_diff(p = p) %>%
#               data.frame(
#                 p = p,
#                 X = X,
#                 N = N,
#                 method = method,
#                 rep = r,
#                 diff = .
#               )
#           }
#         )
#       }
#     )
#   }
# )
#
# diff.df %>%
#   subset(rep < 20) %>%
#   ggplot(aes(p,diff,colour = method,fill = method))+
#   facet_grid(X~method)+
#   geom_jitter(aes(group = interaction(method,p)),alpha = 0.1)+
#   geom_boxplot(aes(group = interaction(method,p)),alpha = 0.2)+
#   theme_bw()
#
#
#
#
# set.seed(42)
# n <- 5L
# node_names <- as.character(1:n)
# # P.ori <- draw_P.ori(n,rng.vec = rnorm(n,50,10),node_names = node_names,min = 0,max = .001,noise.sd = 0.2,adjust.min = 0)
# P.ori <- matrix(0,ncol = n,nrow = n,dimnames = list(node_names,node_names)) %>% {diag(.) <- 0; non.diagonal(.) <- seq(0,1,length.out = n * (n - 1));.}
#
# # distribution from simulation
# ## one single draw
# A.obs <- rbinom_p(N,P.ori)
# A.obs
#
# P.sim <- draw_P.simu(A.obs,N)
# P.sim
#
# A.sim <- rbinom_p(N,P.sim)
# A.sim
#
#
# n.rep <- 200L
# N <- 10L
# A.obs <- rbinom_p(N,P.ori)
# scan.list.test <- replicate(n.rep,simulate_scan(A.obs,N,1:N) %>% get_scanList("A") %>% Reduce("+",.))
# bbinom.test <- replicate(n.rep,simulate_bbinom(A.obs,N))
#
# bbinom.vs.scanlist.df <-
#   rbind_lapply(
#     1:n.rep,
#     function(r) {
#       rbind(
#         expand.grid(i = 1:n,j = 1:n) %>% {
#           .$rep <- r
#           .$N <- N
#           .$method <- factor("scan.list",levels = c("scan.list","bbinom"))
#           .$a <- scan.list.test[,,r][cbind(.$i,.$j)]
#           .
#         },
#         expand.grid(i = 1:n,j = 1:n) %>% {
#           .$rep <- r
#           .$N <- N
#           .$method <- factor("bbinom",levels = c("scan.list","bbinom"))
#           .$a <- bbinom.test[,,r][cbind(.$i,.$j)]
#           .
#         }
#       )
#     }
#   )
#
# bbinom.vs.scanlist.df %>%
#   ggplot(aes(method,a / N,colour = method,fill = method))+
#   facet_grid(i~j)+
#   geom_hline(yintercept = c(0,1),colour = "black")+
#   geom_jitter(alpha = 0.1)+
#   geom_boxplot(alpha = 0.3)+
#   cowplot::theme_half_open()+scale_y_continuous(expand = expansion(mult = c(0,0.1)))
#
# ## distribution with many draws
#
# obs.vs.sim.df <-
#   rbind_lapply(
#     seq(10,200,by = 10),
#     function(N) {
#       A.obs <- rbinom_p(N,P.ori)
#
#       rbind_lapply(
#         1:100,
#         function(r) {
#           Simu <- simulate_scan(A.obs,N)
#
#           # P.sim <- Simu %>% get_scanList("P.obs") %>% median
#
#           A.sim <- Simu %>% get_scanList("A") %>% Reduce("+",.)
#
#           rbind(
#             expand.grid(i = 1:n,j = 1:n) %>% {
#               .$rep <- r
#               .$p <- P.ori[cbind(.$i,.$j)]
#               .$N <- N
#               .$stage <- factor("observed",levels = c("observed","simulated"))
#               .$a <- A.obs[cbind(.$i,.$j)]
#               .
#             },
#             expand.grid(i = 1:n,j = 1:n) %>% {
#               .$rep <- r
#               .$p <- P.ori[cbind(.$i,.$j)]
#               .$N <- N
#               .$stage <- factor("simulated",levels = c("observed","simulated"))
#               .$a <- A.sim[cbind(.$i,.$j)]
#               .
#             }
#           )
#         }
#       )
#     }
#   )
#
# # obs.vs.sim.df %>%
# #   subset(i != j) %>%
# #   ggplot(aes(N,a / N,colour = stage, fill = stage))+
# #   facet_grid(i~j)+
# #   geom_hline(yintercept = c(0,1),colour = "black")+
# #   geom_jitter(aes(group = interaction(stage,N)),alpha = 0.1,height = 0)+
# #   geom_boxplot(aes(group = interaction(stage,N)),alpha = 0.3)+
# #   geom_smooth()+
# #   geom_hline(data = subset(obs.vs.sim.df,stage == "observed"),aes(yintercept = p,colour = stage),lty = "dashed")+
# #   cowplot::theme_half_open()+scale_y_continuous(expand = expansion(mult = c(0,0.1)))
#
# obs.vs.sim.df %>%
#   setDT %>% {
#     .[i != j,by = .(i,j,p,N,stage),
#       .(
#         a = median(a),
#         lower = quantile(a,0.025),
#         upper = quantile(a,0.975)
#       )]
#   } %>%
#   # subset(stage == "observed") %>%
#   ggplot(aes(N,a / N,colour = stage, fill = stage))+
#   facet_grid(i~j)+
#   geom_hline(yintercept = c(0,1),colour = "black")+
#   geom_line(data = subset(obs.vs.sim.df, stage == "observed"),lty = "dotted",alpha = 0.8)+
#   geom_point(data = subset(obs.vs.sim.df, stage == "observed"),alpha = 0.8)+
#   geom_ribbon(aes(ymin = lower / N,ymax = upper / N),alpha = 0.1,colour = NA)+
#   geom_hline(aes(yintercept = p,colour = stage),lty = "dashed")+
#   cowplot::theme_half_open()+scale_y_continuous(expand = expansion(mult = c(0,0.1)))
#
# ## distribution with bootstrap
#
# obs.vs.sim.boot.df <-
#   rbind_lapply(
#     seq(10,100,by = 20),
#     function(N) {
#       A.obs <- rbinom_p(N,P.ori)
#       Simu.to.boot <- simulate_scan(A.obs,N)
#       cat("\nN = ",N,"\n")
#       rbind_lapply(
#         1:500,
#         function(r) {
#           cat("\nreplicate = ",r)
#           A.to.boot <- Simu.to.boot %>% get_scanList("A") %>% Reduce("+",.)
#           A.boot <- Simu.to.boot %>% get_scanList("A") %>% sample(replace = TRUE) %>% Reduce("+",.)
#
#           Simu <- simulate_scan(A.obs,N) %>% get_scanList("A")
#
#           A.sim <- Simu %>% Reduce("+",.)
#           cat("   Starting bootstrap")
#           # A.both <- Simu %>% sample(x = .,replace = TRUE) %>% Reduce("+",.)
#           A.both <- Simu %>% {replicate(50,sample(x = .,replace = TRUE) %>% Reduce("+",.),simplify = FALSE)} %>% Reduce("+",.) / 50
#           cat("\rBootstrap done!            \n")
#
#           rbind(
#             expand.grid(i = 1:n,j = 1:n) %>% {
#               .$rep <- r
#               .$p <- P.ori[cbind(.$i,.$j)]
#               .$N <- N
#               .$stage <- factor("observed",levels = c("observed","simulated","to.boot","bootstrapped","both"))
#               .$a <- A.obs[cbind(.$i,.$j)]
#               .
#             },
#             expand.grid(i = 1:n,j = 1:n) %>% {
#               .$rep <- r
#               .$p <- P.ori[cbind(.$i,.$j)]
#               .$N <- N
#               .$stage <- factor("to.boot",levels = c("observed","simulated","to.boot","bootstrapped","both"))
#               .$a <- A.to.boot[cbind(.$i,.$j)]
#               .
#             },
#             expand.grid(i = 1:n,j = 1:n) %>% {
#               .$rep <- r
#               .$p <- P.ori[cbind(.$i,.$j)]
#               .$N <- N
#               .$stage <- factor("simulated",levels = c("observed","simulated","to.boot","bootstrapped","both"))
#               .$a <- A.sim[cbind(.$i,.$j)]
#               .
#             },
#             expand.grid(i = 1:n,j = 1:n) %>% {
#               .$rep <- r
#               .$p <- P.ori[cbind(.$i,.$j)]
#               .$N <- N
#               .$stage <- factor("bootstrapped",levels = c("observed","simulated","to.boot","bootstrapped","both"))
#               .$a <- A.boot[cbind(.$i,.$j)]
#               .
#             },
#             expand.grid(i = 1:n,j = 1:n) %>% {
#               .$rep <- r
#               .$p <- P.ori[cbind(.$i,.$j)]
#               .$N <- N
#               .$stage <- factor("both",levels = c("observed","simulated","to.boot","bootstrapped","both"))
#               .$a <- A.both[cbind(.$i,.$j)]
#               .
#             }
#           )
#         }
#       )
#     }
#   )
#
#
# obs.vs.sim.boot.df %>%
#   setDT %>% {
#     .[i != j,by = .(i,j,p,N,stage),
#       .(
#         a = median(a),
#         lower = quantile(a,0.025),
#         upper = quantile(a,0.975)
#       )]
#   } %>%
#   ggplot(aes(N,a / N,colour = stage, fill = stage))+
#   facet_grid(i~j)+
#   geom_hline(yintercept = c(0,1),colour = "black")+
#   geom_linerange(aes(ymin = lower / N,ymax = upper / N),alpha = 0.3,position = position_dodge(5))+
#   geom_line(alpha = 1,position = position_dodge(5),lty = "dotted")+
#   geom_point(alpha = 1,position = position_dodge(5))+
#   geom_hline(data = subset(obs.vs.sim.df,stage %in% c("observed","to.boot")),aes(yintercept = p,colour = stage),lty = "dashed")+
#   cowplot::theme_half_open()+scale_y_continuous(expand = expansion(mult = c(0,0.1)))
#
#
# # draw_P.ori <- function(n,p_ij = rnorm(n,mean = 0.5,sd = 0.2),node_names) {
# #   p_ij <- ifelse(p_ij < 0,0,ifelse(p_ij > 1,1,p_ij))
# #   P.ori <- matrix(
# #     p_ij,nrow = n,ncol = n,
# #     dimnames = list(node_names,node_names)
# #   )
# #   diag(P.ori) <- 0
# #   P.ori
# # }
# #
# # binom.mat <- function(n,p.ij = 0.5,dimnames = NULL) {
# #   matrix(rbinom(n * n, 1, prob = p.ij), ncol = n, nrow = n,dimnames = dimnames)
# # }
# #
# # draw_NEff.scans <- function(NEff,n,P.ori,node_names) {
# #   lapply(
# #     1:NEff,
# #     function(i) {
# #       binom.mat(n,P.ori,dimnames = list(node_names,node_names))
# #     }
# #   )
# # }
# #
# # generate_Orig <- function(NEff,n,P.ori,node_names) {
# #   scan.list <- draw_NEff.scans(NEff,n,P.ori,node_names)
# #   Reduce("+",scan.list)
# # }
# #
# # generate_P.approx <- function(NEff,n,P.ori,node_names,method = c("SimuNet","simple")) {
# #   Orig <- generate_Orig(NEff,n,P.ori,node_names)
# #   P.approx <- list(
# #     SimuNet = Orig * (1 - 2/NEff)/NEff + 1/NEff,
# #     simple = Orig / NEff
# #   )
# #   diag(P.approx$SimuNet) <- diag(P.approx$simple) <- 0
# #   P.approx
# # }
# #
# # generate_cor.int.df <- function(NEff,n,P.ori,node_names) {
# #   rbind_lapply(
# #     seq_along(NEff),
# #     function(l) {
# #       neff <- NEff[l]
# #       P.approx <- generate_P.approx(neff,n,P.ori,node_names)
# #       ct.SimuNet <- cor.test(
# #         non.diagonal(P.ori,output = "vector"),
# #         non.diagonal(P.approx$SimuNet,output = "vector")
# #       )
# #       ct.simple <- cor.test(
# #         non.diagonal(P.ori,output = "vector"),
# #         non.diagonal(P.approx$simple,output = "vector")
# #       )
# #       cat("\rCalculating for NEff = ",neff," (",l,"/",length(seq_along(NEff)),")",sep = "")
# #       data.frame(
# #         n = n,
# #         NEff = neff,
# #         method = factor(c("SimuNet","simple")),
# #         inf = c(ct.SimuNet$conf.int[1],ct.simple$conf.int[1]),
# #         sup = c(ct.SimuNet$conf.int[2],ct.simple$conf.int[2])
# #       )
# #     }
# #   )
# # }
#
# # set.seed(42)
# # n <- 15L
# # NEff <- 5L
# # node_names <- as.character(1:n)
# #
# #
# #
# #
# # cor.int <- rbind_lapply(
# #   c(5,10,20,50),
# #   function(n) {
# #     cat("Starting with n =",n,"\n")
# #     node_names <- as.character(1:n)
# #     P.ori <- draw_P.ori(n,runif(n),node_names)
# #     c.int <-
# #       rbind_lapply(
# #       1:10,
# #       function(r) {
# #         data.frame(
# #          r = r,
# #          generate_cor.int.df(c(2:10 %o% 10^(1:3)),n,P.ori,node_names)
# #         )
# #       }
# #     )
# #     cat("\n")
# #     c.int
# #   }
# # )
# #
# # cor.int <- data.table(cor.int)
# # summary(cor.int)
# #
# # cor.int[,cor := (inf + sup) / 2,by = .(n,NEff)]
# # cor.summ <- cor.int[,.(cor = mean(cor)),by = .(n,NEff,method)]
# # cor.summ$inf <- cor.int[cor.int[, .I[inf == min(inf)], by = .(n,NEff,method)]$V1]$inf
# # cor.summ$sup <- cor.int[cor.int[, .I[sup == max(sup)], by = .(n,NEff,method)]$V1]$sup
# #
# # ggplot(cor.summ,aes(NEff,cor,colour = method,group = method))+
# #   facet_grid(.~n)+
# #   geom_linerange(aes(ymin = inf,ymax = sup),alpha = 0.4)+#,colour = "tomato")+
# #   geom_line(lty="dashed",alpha = .7)+
# #   geom_point()+#colour = "tomato")+
# #   scale_x_log10(breaks = scales::log_breaks(n = 20))+
# #   # scale_color_discrete()+
# #   theme_bw()
#
#
#
# # draw_greg.index <- function(n,rng.vec = rnbinom(n,mu = 2,size = .7),min = 0.025,max = 0.975,node_names = NULL) {
# #   names(rng.vec) <- node_names
# #   scales::rescale(rng.vec,to = c(min,max))
# # }
# #
# # shape_dyads.df <- function(n,greg.index) {
# #   dyads <- expand.grid(i = 1:n,j = 1:n)
# #   greg.prod <- greg.index %o% greg.index
# #   diag(greg.prod) <- 0
# #   dyads$greg.prod <- greg.prod[cbind(dyads$i,dyads$j)]
# #   dyads
# # }
# #
# # draw_dyadic.noise <- function(dyads.df,adjust.min = 0.025,mean = 0, sd = 0.25) {
# #   truncnorm::rtruncnorm(n = nrow(dyads.df),
# #                         a = -dyads.df$greg.prod + adjust.min,
# #                         b = 1 - dyads.df$greg.prod - adjust.min,
# #                         mean = mean,
# #                         sd = sd
# #   )
# # }
# #
# # add_dyadic.noise <- function(dyads.df,adjust.min = 0.025,mean = 0, sd = 0.25) {
# #   dyadic.noise <- draw_dyadic.noise(dyads.df,adjust.min,mean,sd)
# #   dyads.df$dyadic.noise <- ifelse(dyads.df$i != dyads.df$j,dyadic.noise,0)
# #   dyads.df
# # }
# #
# # draw_P.ori <- function(n,node_names,greg.index.args,dyadic.noise.args) {
# #   # Randomize an individual sociality attribute
# #   greg.index <- draw_greg.index(n,rng.vec,
# #                                 min = min,max = max,
# #                                 node_names = node_names
# #   )
# #
# #   # Randomize specific dyad affinities
# #   ## First simple deterministic product of greg.index
# #   dyads.df <- shape_dyads.df(n,greg.index)
# #
# #   ## add a uniform noise around the product
# #   dyads.df <- add_dyadic.noise(dyads.df,adjust.min = adjust.min,mean = mean, sd = sd)
# #   dyads.df$p.ij <- dyads.df$greg.prod + dyads.df$dyadic.noise
# #
# #   # generate P.ori
# #   matrix(dyads.df$p.ij,nrow = n,ncol = n,dimnames = list(node_names,node_names))
# # }
# #
#
#
# # ## Parameters
# # n <- 15L
# # node_names <- letters[1:n]
# #
# # greg.index.args <- list(
# #   rng.vec = rnbinom(n,mu = 2,size = 0.8),
# #   min = 0,
# #   max = 0.8
# # )
# #
# # dyadic.noise.args <- list(
# #   adjust.min = 0,
# #   mean = 0,
# #   sd = 0.001
# # )
# #
# # P.ori <- draw_P.ori(n,node_names,greg.index.args,dyadic.noise.args)
# # P.ori
# #
# #
# #
# # P.ori.list <- replicate(5,draw_P.ori(n,node_names,greg.index.args,dyadic.noise.args),simplify = FALSE)
# # P.ori.multiarray <- abind::abind(P.ori.list,along = 3)
# # object.size(P.ori.list)
# # object.size(P.ori.multiarray)
# #
# # bench <- microbenchmark::microbenchmark(
# #   lapply = lapply(
# #     P.ori.list,
# #     function(m) {
# #       m[6:8,1:5] + m[4:6,6:10]
# #     }
# #   ),
# #   abind = {P.ori.multiarray[6:8,1:5,1:5] + P.ori.multiarray[4:6,6:10,1:5]},
# #   times = 50,unit = "ns"
# # )
# #
# # boxplot(bench)
# #
# # test <- lapply(
# #   P.ori.list,
# #   function(m) {
# #     m[6:8,1:5] + m[4:6,6:10]
# #   }
# # )
# # tost <- P.ori.multiarray[6:8,1:5,1:5] + P.ori.multiarray[4:6,6:10,1:5]
# #
# # sapply(
# #   seq_along(test),
# #   function(i) {
# #     identical(test[[i]],tost[,,i])
# #   }
# # )
#
#
#
# # set.seed(42)
# # nE <- 10L
# # P.test <- runif(1)
# # P.test
# # O.test <- rbinom(1,nE,P.test)
# # O.test
# # P.hat <- O.test / nE
# # P.hat
# #
# #
# # set.seed(42)
# # nE <- 10L
# # P.test <- runif(1)
# # P.test
# # O.test <- rbinom(1,nE,P.test)
# # O.test
# # P.hat <- O.test / nE
# # P.hat
# # CI <- P.hat + c(-1.96,1.96) * sqrt((P.hat * (1 - P.hat)) / nE)
# # CI
# #
# # binom::binom.confint(O.test,nE,conf.level = 0.95)
# #
# # determine_CI <- function(p,n,alpha = 0.05) {
# #   o <- rbinom(1,n,p)
# #   p.hat <-  o / n
# #   p.hat %>%
# #     {. + c(-1,1) * qnorm(1 - (alpha / 2)) * sqrt((. * (1 - .)) / n)} %>%
# #     data.frame(p = p,n = n,o = o,p.hat = p.hat,lower = .[1],upper = .[2])
# # }
# #
# # CI.df <- ConfiNet::rbind_lapply(
# #   c(1:50 * 10),
# #   function(n) {
# #     binom::binom.confint(rbinom(1,n,P.test),n,conf.level = 0.95)
# #   }
# # )
# #
# # library(ggplot2)
# # CI.df %>%
# #   subset(method %in% c("agresti-coull","bayes","exact","wilson")) %>%
# #
# # ggplot(aes(n,mean,fill = method))+
# #   facet_grid(method~.)+
# #   geom_ribbon(aes(ymin = lower,ymax = upper),alpha = 0.2)+
# #   geom_line()+
# #   geom_point(shape = 21,fill = "white")+
# #   geom_hline(yintercept = P.test,lty = "dashed",colour = "darkred",alpha = .6)+
# #   theme_bw()
# #
# # set.seed(42)
# # nE <- 20L
# # # np <- 5L
# # # P.test <- sort(runif(np))
# # P.test <- seq(from = 0,to = 1,by = 0.20)
# # np <- length(P.test)
# #
# # rbind_lapply(
# #   1:50 * 10,
# #   function(n) {
# #     rbind_lapply(
# #       1:50,
# #       function(r) {
# #         P.test %>%
# #           rbinom(n = length(.),size = n,prob = .) %>%
# #           binom::binom.confint(n = n,methods = c("wilson","exact")) %>%
# #           cbind(P = P.test,Probability = paste0("p = ",round(P.test,2)),replicate = r,.)
# #       }
# #     )
# #   }
# # ) %>%
# #   ggplot(aes(n,mean,fill = Probability,colour = Probability,group = interaction(replicate,Probability)))+
# #   facet_grid(method~.)+
# #   geom_ribbon(aes(ymin = lower,ymax = upper),alpha = 0.01,colour = NA)+
# #   geom_line(alpha = 0.04)+
# #   geom_point(shape = 21,alpha = 0.02)+
# #   geom_hline(aes(yintercept = P,colour = Probability),lty = "dashed",size = 2)+
# #   geom_hline(yintercept = 0,colour = "black")+
# #   geom_vline(xintercept = 10,colour = "black")+
# #   cowplot::theme_half_open()+scale_x_continuous(expand = expansion(mult = c(0,0.01)))+scale_y_continuous(expand = expansion(mult = c(0,0.1)))
# #
# # n <- 1000L
# # alpha.df <-
# #   rbind_lapply(
# #   1:999 / 1000,
# #   function(alpha) {
# #     rbind_lapply(
# #       1:10,
# #       function(r) {
# #         P.test[2] %>%
# #           rbinom(n = length(.),size = n,prob = .) %>%
# #           binom::binom.confint(n = n,methods = c("wilson","exact"),conf.level = 1 - alpha) %>%
# #           cbind(P = P.test,Probability = paste0("p = ",round(P.test[2],2)),alpha = alpha,replicate = r,.)
# #       }
# #     )
# #   }
# # )
# # alpha.df %>%
# #   ggplot(aes(mean))+
# #   facet_grid(method~.)+
# #   geom_histogram(binwidth = 0.002,fill = "grey80",colour = "grey30", alpha = 0.8)+cowplot::theme_half_open()#+
# #
# #   # ggplot(aes(alpha,mean,fill = Probability,colour = Probability,group = interaction(replicate,Probability)))+
# #   # facet_grid(method~.)+
# #   # geom_ribbon(aes(ymin = lower,ymax = upper),alpha = 0.01,colour = NA)+
# #   # geom_line(alpha = 0.04)+
# #   # geom_point(shape = 21,alpha = 0.02)+
# #   # cowplot::theme_half_open()#+
# #   # scale_x_continuous(expand = expansion(mult = c(0,0.01)))+scale_y_continuous(expand = expansion(mult = c(0,0.1)))
# #
# # library(dplyr)
# # library(ggplot2)
# #
# # calculate_CIs <- function(p,N,
# #                           alphas = c(
# #                             seq(from = 0.01,
# #                                        to = 0.99,
# #                                        by = 0.01
# #                             ),
# #                             seq(.990,.999,by = 0.001)
# #                           )
# #                           ) {
# #   X <-
# #     p %>%
# #     rbinom(n = length(p),size = N,prob = .)
# #
# #   do.call(
# #     rbind,
# #     lapply(
# #       alphas,
# #       function(alpha) {
# #         X %>%
# #           binom::binom.confint(n = N,methods ="exact",conf.level = 1 - alpha) %>%
# #           cbind(p = p,alpha = alpha,.)
# #       }
# #     )
# #   )
# # }
# #
# # set.seed(1)
# # p1 <- seq(from = .15,to = 1,by = 0.1)
# #
# # df <- rbind(calculate_CIs(p = p1,N = 10),
# #       calculate_CIs(p = p1,N = 50),
# #       calculate_CIs(p = p1,N = 100)) %>%
# #   mutate(N = factor(
# #     paste0("N = ",n),
# #     levels = paste0("N = ",c("10","50","100"))
# #   )
# #   )
# #
# # df %>%
# #   subset(n == 10 & p == 0.15) %>%
# #   ggplot(aes(mean,alpha,colour = N))+
# #   facet_grid(N ~  paste0("p1 = ",p))+
# #   geom_errorbarh(aes(xmin = lower,xmax = upper))+
# #   geom_vline(aes(xintercept = mean),lty = "dotted",size = 1.1)+
# #   geom_vline(aes(xintercept = p),colour = "red",lty = "dashed",size = 1.1)+
# #   theme_bw()+
# #   xlab("p1_hat")+
# #   guides(colour = FALSE)+
# #   scale_x_continuous(breaks = c(0,0.5,1),limits = c(0,1))+
# #   scale_y_continuous(expand = expansion(mult = c(0,0.01)),limits = c(0,1))+
# #   geom_function(fun = function(p) (p^1 * (1 - p)^(10-1)) / dbinom(1,10,0.1))
# #
# #
# # posterior_fun <- function(p,X,N) {
# #   (p^X * (1 - p)^(N - X)) / pbinom(X,N,X/N)
# # }
# #
# # integrate(function(p) posterior_fun(p,1,10),lower = 0,upper = 1,)
# #
# # ggplot(df,aes(mean))+
# #   geom_function(fun = function(mean) posterior_fun(mean,1,10))+
# #   theme_bw()+
# #   guides(colour = FALSE)+
# #   scale_x_continuous(breaks = c(0,0.5,1),limits = c(0,1))+
# #   scale_y_continuous(expand = expansion(mult = c(0,0.01)))
# # (test <- runif(1))
# # test %>%
# #   posterior_fun(1,10)
# #
# # N <- 10L
# # p <- seq(0,1,by = 0.01)
# # X <- rbinom(1,N,0.3)
# # plot(p,p^(X) * (1 - p)^(N - X),type = "l")
# #
# #
# # p <- seq(0,1, by=.01)
# # like = p^60*(1-p)^40
# # mle = p[like==max(like)];  mle
# # plot(p, like, type="l")
# # abline(h=0, col="green2")
# # abline(v = mle, col="red", lwd=2)
# #
# #
# # df <-
# #   do.call(
# #     rbind,
# #     lapply(
# #       seq(from = 0.01,to = 0.99,by = 0.01),
# #       function(alpha) {
# #         do.call(
# #           rbind,
# #           lapply(
# #             1:1,
# #             function(r) {
# #               p %>%
# #                 rbinom(n = length(p),size = n,prob = .) %>%
# #                 binom::binom.confint(
# #                   n = n,
# #                   methods = c("agresti-coull",
# #                               "wilson",
# #                               "exact"
# #                   ),
# #                   conf.level = 1 - alpha
# #                 ) %>%
# #                 cbind(
# #                   p = p,
# #                   alpha = alpha,
# #                   replicate = r,
# #                   .
# #                 )
# #             }
# #           )
# #         )
# #       }
# #     )
# #   )
# # df %>%
# #   ggplot(aes(mean,fill = method,colour = method))+
# #   facet_grid(. ~ p + method)+
# #   geom_histogram(bins = 30,alpha  = .6)+
# #   theme_bw()
# # library(data.table)
# # setDT(df)
# # df[,.(mean = mean(mean),se.mean = sd(mean)),by = .(p,method,replicate)]
# #
#
#
# #
# #
# # library(ggplot2)
# # set.seed(42)
# # param <-
# #   expand.grid(
# #     p = seq(0,1,.2),
# #     N = c(10,20,50,100)
# #   )
# #
# # param$group <- as.character(1:nrow(param))
# #
# # df <- rbind_lapply(
# #   1:nrow(param),
# #   function(i) {
# #     p <- param$p[i]
# #     N <- param$N[i]
# #     group <- param$group[i]
# #
# #     X <- rbinom(1,N,p)
# #
# #     p.hat <- rbeta(500,X + 1, (N - X) + 1)
# #
# #     data.table(
# #       p = p,
# #       N = N,
# #       X = X,
# #       shape1 = X + 1,
# #       shape2 = N - X + 1,
# #       group = group,
# #       p.hat = p.hat
# #     )
# #   }
# # )
# #
# # df.summary <- df[,.(p.hat = mean(p.hat)),by = .(p,N,X,group,shape1,shape2)]
# #
# # df %>%
# #   ggplot(aes(p.hat,fill = as.factor(p),colour = as.factor(p),group = group))+
# #   facet_grid(N~p,scales = "free_y")+
# #   geom_histogram(aes(y = stat(density)),colour = NA,binwidth = 0.01,alpha = 0.4)+
# #   lapply(
# #     1:nrow(df.summary),
# #     function(i) {
# #       geom_function(data = df.summary[i,],
# #                     aes(group = group),
# #                     lty = "dashed",
# #                     size = .8,
# #                     # colour = "grey60",
# #                     fun = function(x) dbeta(x,
# #                                             shape1 = df.summary$shape1[i],
# #                                             shape2 = df.summary$shape2[i]
# #                     )
# #       )
# #     }
# #   )+
# #   geom_vline(aes(xintercept = p),colour = "red")+
# #   geom_text(data = df.summary,aes(0.5,10,label = paste0("X = ",X,"\np.hat = ",round(p.hat,2),"\nB(",shape1,",",shape2,")")))+
# #   geom_vline(aes(xintercept = 0),colour = "black")+
# #   geom_hline(aes(yintercept = 0),colour = "black")+
# #   geom_vline(aes(xintercept = X / N),colour = "green")+
# #   geom_vline(data = df.summary,aes(xintercept = p.hat),colour = "blue")+
# #   scale_x_continuous(breaks = c(0,0.5,1),labels = c("0","0.5","1"),limits = c(0,1),expand = expansion(mult = c(0,0.01)))+
# #   scale_y_continuous(expand = expansion(mult = c(0,0.01)))+
# #   cowplot::theme_half_open()
# #
# #
# # ggplot(aes(p.hat,fill = as.factor(p),colour = as.factor(p)))+
# #   facet_grid(N~p)+
# #   geom_histogram(binwidth = 0.01,alpha = 0.7)+
# #   geom_vline(aes(xintercept = p))+
# #   geom_function(data = subset(df,p = 0.2),fun = function(p) dbeta(p,X + 1, N - X + 1))+
# #   cowplot::theme_half_open()
#
#

# source library and custom functions -------------------------------------
source(".WIP/validation_tools.R")

set.seed(42)
n <- 5L
N <- 10L
n.rep <- 5L

P <- generate_P.seq(n)
P.df <- extract_xij(P,x.name = "p")

P %>%
  draw_scanList(n.scans = N) %>%
  sum_scanList %>%
  extract_xij(N = N,x.name = "a",rep = r) %>%
  merge(P.df,.,by = c("i","j"))


test.sl <-
  P %>%
  draw_scanList(n.scans = N)

test.sl %>%
  bootstrap_scanList(n.rep) %>%
  calculate_SNm(SNm_fun = function(mat) DirectedClustering::ClustBCG(mat,type = "directed")$GlobaltotalCC,Xapply = sapply) %>%
  data.table(r = 1:n.rep,CC = .) %>%
  ggplot(aes("",CC))+
  geom_jitter(alpha = 0.2,colour = "royalblue")+
  geom_boxplot(alpha = 0.5,colour = "royalblue",fill = "royalblue")+
  theme_cowplot()

calculate_EV <- function(graph) {igraph::eigen_centrality(graph,directed = TRUE)$vector}
calculate_PR <- function(graph) {igraph::page.rank(graph,directed = TRUE,weights = igraph::E(graph)$weight)$vector}

test.sl %>%
  bootstrap_scanList(n.rep) %>%
  calculate_SNm(SNm_fun = calculate_PR,mode = "directed",transform.to.igraph = TRUE,Xapply = function(X,FUN,...) do.call(rbind,lapply(X,FUN,...))) %>%
  data.table %>%
  cbind(r = 1:n.rep,.) %>%
  melt.data.table(id.vars = 1,measure.vars = 1:n + 1,variable.name = "node",value.name = "PR") %>%
  ggplot(aes(node,PR,colour = node, fill = node))+
  geom_jitter(alpha = 0.2)+
  geom_boxplot(alpha = 0.5)+
  theme_cowplot()

test.sl %>%
  bootstrap_scanList(n.rep) %>%
  calculate_SNm(SNm_fun = calculate_PR,mode = "directed",transform.to.igraph = TRUE,Xapply = function(X,FUN,...) do.call(rbind,lapply(X,FUN,...))) %>%
  data.table %>%
  cbind(r = 1:n.rep,.) %>%
  melt.data.table(id.vars = 1,measure.vars = 1:n + 1,variable.name = "node",value.name = "PR") %>% {
    .[,by = .(node),.(PR = median(PR),lower = quantile(PR,probs = 0.025),upper = quantile(PR,probs = 0.975),n.rep = .N)]
  }

P %>%
  igraph::graph.adjacency(mode = "directed",weighted = TRUE) %>%
  calculate_PR


lapply(
  c(10,50,100,1000),
  function(N) {
    lapply(
      1:30,
      function(r) {
        P %>%
          draw_scanList(n.scans = N) %>%
          sum_scanList %>%
          extract_xij(N = N,x.name = "a",rep = r) %>%
          merge(P.df,.,by = c("i","j"))
      }
    ) %>% do.call(rbind,.)
  }
) %>%
  do.call(rbind,.) %>%
  ggplot(aes(p,a / N))+
  facet_grid(N~.)+
  geom_point(alpha = 0.1)+
  theme_cowplot()


P %>%
  mat_rbinom(N) %>% {
    rbind(
      replicate(
        n.rep,
        simplify = FALSE,
        mat_rbbinom(.,N) %>%
          extract_xij(N = N,x.name = "a",method = "Beta-binomial") %>%
          merge(P.df,.,by = c("i","j")) %>% data.table
      ),
      replicate(
        n.rep,
        simplify = FALSE,
        draw_scanList(A = .,N = N) %>%
          sum_scanList() %>%
          extract_xij(N = N,x.name = "a",method = "SimuNet") %>%
          merge(P.df,.,by = c("i","j")) %>% data.table
      )
    )
  } %>%
  do.call(rbind,.) %>%
  ggplot(aes(p,a / N,colour = method,fill = method))+
  geom_jitter(alpha = 0.03)+
  geom_boxplot(aes(group = interaction(method,p)),alpha = 0.5)+
  geom_smooth(method = "lm")+
  theme_cowplot()

bbinom.boot_across.N.df <-
  lapply(
    c(10,20,50,100),
    function(N) {
      lapply(
        1:100,
        function(r) {
          scan.list <- P %>% draw_scanList(N = N)
          rbind(
            replicate(
              n.rep,
              simplify = FALSE,
              scan.list %>%
                sum_scanList %>%
                mat_rbbinom(.,N) %>%
                extract_xij(N = N,x.name = "a",r = r,method = "Beta-binomial") %>%
                merge(P.df,.,by = c("i","j")) %>% data.table
            ) %>% do.call(rbind,.),
            replicate(
              n.rep,
              simplify = FALSE,
              scan.list %>%
                resample_scanList() %>%
                sum_scanList() %>%
                extract_xij(N = N,x.name = "a",r = r,method = "Bootstrap") %>%
                merge(P.df,.,by = c("i","j")) %>% data.table
            ) %>% do.call(rbind,.)
          )
        }
      ) %>% do.call(rbind,.)
    }
  ) %>% do.call(rbind,.)

bbinom.boot_across.N.df %>%
  ggplot(aes(p,a / N,colour = method,fill = method))+
  facet_grid(N~.)+
  # geom_jitter(alpha = 0.03)+
  geom_boxplot(aes(group = interaction(method,p)),alpha = 0.5)+
  # geom_smooth(method = "lm")+
  theme_cowplot()


bbinom.boot_across.N.df[,abs.err := abs(p - a / N),]
bbinom.boot_across.N.df %>% {
  .[,by = .(i,j,p,N,method),.(abs.err = median(abs.err),lower = quantile(abs.err,0.025),upper = quantile(abs.err,0.975))]
} %>%
  ggplot(aes(N,abs.err,fill = method,color = method))+
  facet_grid(i~j)+
  geom_ribbon(aes(ymin = lower, ymax = upper),colour = NA,alpha = 0.2)+
  geom_line(alpha = 1)+
  theme_cowplot()

bbinom.boot_across.N.df[,by = .(N,i,j,p,method,r),c("p.lower","p.upper") := .(quantile(a / N,0.025),quantile(a / N,0.975))]

bbinom.boot_across.N.df[,c("is.lower","is.upper") := .(as.integer(p < p.lower),as.integer(p > p.upper)),]

bbinom.boot_across.N.df[is.lower == 1 & method == "Beta-binomial",.N,by = .(i,j,p,r,method,N,a)]

bbinom.boot_across.N.df[is.lower == 1 & method == "Beta-binomial",by = .(i,j,p,r,method,N),.N][,by = .(i,j,p,method,N),.N / N]


bbinom.boot_across.N.df[,by = .(N,i,j,p,method,r),.(false.pos = sum(.SD$is.lower)/.N,false.neg = sum(.SD$is.upper)/.N)] %>% {
  .[,by = .(N,i,j,p,method),
    .(
      false.neg = median(false.neg),
      lower = quantile(false.neg,0.025),
      upper = quantile(false.neg,.975)
    )
  ]
} %>%
  ggplot(aes(N,false.neg,fill = method,colour = method))+
  facet_grid(i~j)+
  geom_ribbon(aes(ymin = lower, ymax = upper),colour = NA,alpha = 0.2)+
  geom_line()+
  theme_cowplot()


run_experiment <- function(P,N,n.rep = 100,n.distrib.sample = 100) {
  lapply(
    1:n.rep,
    function(r) {
      scan.list <- P %>% draw_scanList(N = N)
      A <- scan.list %>% sum_scanList
      rbind(
        replicate(
          n.distrib.sample,
          simplify = FALSE,
          A %>%
            mat_rbbinom(.,N) %>%
            extract_xij(N = N,x.name = "a",r1 = r,method = "Beta-binomial") %>%
            merge(P.df,.,by = c("i","j")) %>% data.table
        ) %>% do.call(rbind,.),
        replicate(
          n.distrib.sample,
          simplify = FALSE,
          scan.list %>%
            resample_scanList() %>%
            sum_scanList() %>%
            extract_xij(N = N,x.name = "a",r1 = r,method = "Bootstrap") %>%
            merge(P.df,.,by = c("i","j")) %>% data.table
        ) %>% do.call(rbind,.)
      )
    }
  ) %>% do.call(rbind,.)
}

test.dt <- run_experiment(P,N,5,5)

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

# Loading used packages --------------------------------------
library(ggplot2)
library(cowplot)
library(ggridges)
library(data.table)
library(dplyr)
library(extraDistr)
library(binaryLogic)

as.fixed_length.binary <- function(x,length.out,...) {
  if (length(x) == 1) {
    as.binary(x,...) %>%
      fillUpToBit(n = length.out)
  } else {
    as.binary(x,...) %>%
      lapply(.,fillUpToBit,n = length.out)
  }
}

bin2int <- function(x) {
  rev(as.integer(c(x)))
}

# latest ------------------------------------------------------------------
n <- 5L
node_names <- letters[1:n]
M <- matrix(0L,n,n,dimnames = list(node_names,node_names))
N <- 50L
n.non.diag <- (n^2 - n)
n.possib <- 2^n.non.diag
n.possib.und <- 2^(n.non.diag / 2)

p.seq <- seq(0,1,length = n.non.diag / 2)
p.seq <- rbeta(n.non.diag / 2,shape1 = 1,1)
p.seq[3] <- 0.9

# slow for n >= 5
# 1:n.possib %>%
#   {as.fixed_length.binary(. - 1,n^2 - n)} %>%
#   lapply(bin2int) %>%
#   lapply(function(vec) {non.diagonal(M) <- vec;M})

1:n.possib.und %>%
  {as.fixed_length.binary(. - 1,(n^2 - n) / 2)} %>%
  lapply(bin2int) %>%
  lapply(function(vec) {M[lower.tri(M)] <- vec;M}) %>%
  Filter(function(M) !((M[4,1] == 1) & any(M[,2] == 1,M[2,1] == 1)),.) %>%
  Filter(function(M) M[4,1] == 1,.) %>%
  {.[sample(1:length(.),1000,replace = TRUE)]} %>%
  Reduce("+",.)

coords <-
  expand.grid(i = 1:n,j = 1:n) %>%
  data.table %>%
  subset(i > j) %>%
  cbind(.,p = p.seq)

1:n.possib.und %>%
  {as.fixed_length.binary(. - 1,(n^2 - n) / 2)} %>%
  lapply(bin2int) %>%
  lapply(function(vec) cbind(coords,s = vec)) %>%
  lapply(function(df) cbind(df,likelihood = ifelse(df$s == 1,df$p, 1 - df$p))) %>%
  lapply(function(df) list(df = df,likelihood = prod(df$likelihood))) %>%
  Filter(function(df) df$likelihood != 0,.) %>%
  Filter(function(df) !(df$df[i == 4 & j == 1]$s == 1 & any(df$df[j == 2]$s == 1)),.) %>%
  Filter(function(df) !(df$df[i == 4 & j == 1]$s == 1 & df$df[i == 2 & j == 1]$s == 1),.) %>%
  # Filter(function(df) df$df[i == 4 & j == 1]$s == 1,.) %>%
  lapply(function(df) list(M = {M[lower.tri(M)] <- df$df$s;M},likelihood = df$likelihood)) %>%
  {lapply(.,function(df) list(M = df$M,p = df$likelihood / sapply(.,function(df) df$likelihood) %>% sum))} %>%
  {.[sample(1:length(.),1000,replace = TRUE,prob = sapply(.,function(df) df$p))]} %>%
  lapply(function(df) df$M) %>%
  Reduce("+",.) %>%
  igraph::graph.adjacency(mode = "lower",weighted = TRUE) %>%
  {plot.igraph(x = .,edge.width = igraph::edge.attributes(.)$weight/100)}



ceiling(runif(N,0,n.possib)) %>%
  {as.fixed_length.binary(. - 1,n^2 - n)} %>%
  lapply(bin2int) %>%
  lapply(function(vec) {non.diagonal(M) <- vec;M}) %>%
  Reduce("+",.)



non.diagonal(M) <-
  as.binary(42-1) %>% {c(.,rep(0,n^2 - n - length(.)))}
M

# old stuff ---------------------------------------------------------------

set.seed(42)
N <- 50L
X <- seq(0,N,by = 2)
n.rep <- 1000L

test <-
  replicate(
    n = n.rep,
    simplify = FALSE,
    expr = {
      data.table(X = X,
                 method = rep(factor(c("GT","rbbinom","beta.then.binom","scan_sum","bbernouilli","fixed_scan_sum"),levels = c("GT","rbbinom","beta.then.binom","scan_sum","bbernouilli","fixed_scan_sum")),each = length(X)),
                 X.new = c(
                   GT = rbinom(n = length(X),size = N,prob = X / N),
                   rbbinom = rbbinom(n = length(X),size = N,alpha = X + 0.5,beta = N - X + 0.5),
                   beta.then.binom = rbeta(n = length(X),X + 0.5,N - X + 0.5) %>% rbinom(n = length(X),size = N,prob = .),
                   scan_sum = replicate(n = N,rbeta(n = length(X),X + 0.5,N - X + 0.5) %>% rbinom(n = length(X),size = 1,prob = .)) %>% rowSums(),
                   bbernouilli = rbinom(length(X),N,(X + 0.5) / (X + 0.5 + N - X + 0.5)),
                   fixed_scan_sum = rbeta(n = length(X),X + 0.5,N - X + 0.5) %>% {replicate(n = N,rbinom(n = length(X),size = 1,prob = .))} %>% rowSums()
                 )

      )
    }
  ) %>%
    do.call(rbind,.) %>% cbind(rep = rep(1:n.rep,each = length(X)),.)

test %>%
  subset(X %in% seq(0,N,by = 10)) %>%
  ggplot(aes(X.new,interaction(method,X),fill = method))+
  geom_density_ridges()+
  theme_cowplot()

test.summ <-
  test[,.(mean = mean(X.new),median = median(X.new),sd = sd(X.new),lower = quantile(X.new,0.025),upper = quantile(X.new,0.975)),by = .(X,method)][order(X,method)]

test.summ %>%
  ggplot(aes(X,sd,colour = method))+geom_point()+theme_minimal_grid()

test.summ %>%
  ggplot(aes(X,mean,colour = method,fill = method))+geom_point(position = position_dodge(2))+
  geom_linerange(aes(ymin = lower,ymax = upper),position = position_dodge(2))+
  theme_minimal_grid()




# draft -------------------------------------------------------------------
X = 0:50
plot(seq(0,1,0.01),dbeta(seq(0,1,0.01),0.5,0.5),type = "l")

dbbinom(X,50,0.5,0.5)
plot(X,dbbinom(X,50,0.5,0.5))

N <- 50L
data <- rbinom(1,N,0.15)

dbbinom(X,N,data + .5,N - data + .5)
plot(X,dbbinom(X,N,data + .5,N - data + .5))

dbbinom(1,1,data + .5,N - data + .5)

N <- 5000L
data <- rbinom(1,N,0.15)
dbbinom(1,1,data + .5,N - data + .5)


P.to.test <- seq(0,1,.05)
N <- 5000L

get_p.hat_from_bbinom <- function(p,N,n.rep = 100,alpha.prior = 0.5,beta.prior = 0.5) {
  replicate(
    n = n.rep,
    simplify = FALSE,
    N %>%
      rbinom(p,.,prob = p) %>%
      {dbbinom(1,1,. + alpha.prior,N - . + beta.prior)} %>%
      data.table(p.hat = .)
  ) %>%
    do.call(rbind,.) %>%
    cbind(
      p = rep(P.to.test,n.rep),
      N = N,
      rep = rep(1:n.rep,each = length(p)),
      .
    )
}

rbind(
  get_p.hat_from_bbinom(P.to.test,50,alpha.prior = 1,beta.prior = 1),
  get_p.hat_from_bbinom(P.to.test,500,alpha.prior = 1,beta.prior = 1),
  get_p.hat_from_bbinom(P.to.test,5000,alpha.prior = 1,beta.prior = 1)
) %>%
  ggplot(aes(p,p.hat,group = interaction(N,p),colour = factor(N),fill = factor(N)))+
  geom_hline(yintercept = 0)+
  geom_jitter(alpha = 0.05)+
  geom_smooth(inherit.aes = TRUE,aes(group = N),method = "lm")+
  geom_boxplot(alpha = 0.25)+
  theme_cowplot()

set.seed(42)
data <-  50 %>%
  rbinom(1,.,prob = .3)
data
data / 50
data %>% {dbbinom(1,1,. + .5,50 - . + .5)}


50 %>%
  rbinom(1,.,prob = .2) %>% {
    data.table(
      method = factor(rep(c("betabinom","draw.once.beta.prob.then.sum","sum.betabern","redraw.beta.prob.then.sum"),each = 10000)),
      x = c(
        rbbinom(n = 10000,size = 50,alpha = . + .5,beta = 50 - . + .5),
        rbeta(50,. + .5,50 - . + .5) %>% {replicate(n = 10000,rbinom(.,1,.))} %>% colSums,
        replicate(n = 10000,rbbinom(n = 50,size = 1,. + .5,50 - . + .5)  %>% sum),
        replicate(n = 10000,rbeta(50,. + .5,50 - . + .5) %>% {rbinom(.,1,.)} %>% sum)
      )
    )
  } %>%
  ggplot(aes(x,fill = method))+facet_grid(method~.)+
  geom_histogram(binwidth = 1,alpha = 0.7,colour = "black")+theme_cowplot()


rbeta(50,20 + .5,50 - 20 + .5) %>% {replicate(n = 10,rbinom(.,1,.))} %>% colSums


# Loading used packages --------------------------------------
library(ggplot2)
library(cowplot)
library(ggridges)
library(data.table)
library(dplyr)
library(extraDistr)

set.seed(43)
N <- 100L
n.rep <- 1000L
n.obs.draw <- 100L
p.GT <- runif(1)

X.GT <- p.GT %>% rbinom(n.rep,N,.)
dt.GT <- data.table(obs.draw = NA,
                    method = factor(rep(c("GT"),each = n.rep),
                                    levels = c("GT","bootstrap","SimuNet","Beta-Binomial","Sum of Beta-Bernouilli")),
                    rep = 1:n.rep,
                    X.obs = X.GT,
                    X = X.GT
)

full.dt <-
  replicate(
    n = n.obs.draw,simplify = FALSE,
    {
      scan.list.obs <- p.GT %>% {replicate(n = N,rbinom(1,1,.))}
      X.obs <- scan.list.obs %>% sum

      X.bootstrap <- replicate(n = n.rep,scan.list.obs %>% sample(N,TRUE) %>% sum)

      X.SimuNet <- replicate(n = n.rep,rbeta(N,X.obs + 0.5,N - X.obs + 0.5) %>% rbinom(.,1,.) %>% sum)

      X.BBinom <- rbbinom(n.rep,N,X.obs + 0.5,N - X.obs + 0.5)

      X.BBern <- replicate(N,rbbinom(n.rep,1,X.obs + 0.5,N - X.obs + 0.5)) %>% rowSums

      X.BBprob <- replicate(n.rep,dbbinom(1,1,X.obs + 0.5,N - X.obs + 0.5) %>% rbinom(N,1,.) %>% sum)

      data.table(
        method = factor(rep(c("bootstrap","SimuNet","Beta-Binomial","Sum of Beta-Bernouilli","Betabern prob then binom"),each = n.rep)),
        rep = rep(1:n.rep,times = 5),
        X.obs = X.obs,
        X = c(X.bootstrap,X.SimuNet,X.BBinom,X.BBern,X.BBprob)
      )
    }
  ) %>% do.call(rbind,.) %>%
    {cbind(obs.draw = rep(1:n.obs.draw,each = nrow(.) / n.obs.draw),.)} %>%
    rbind(dt.GT,.)


full.dt %>%
  subset(!is.na(obs.draw)) %>%
  subset(obs.draw <= 2) %>%
  ggplot(aes(X / N,colour = method,fill = method))+
  facet_grid(method~.)+
  geom_line(aes(group = obs.draw),stat="density",alpha = 0.15)+
  geom_density(aes(group = obs.draw),colour = NA,alpha = 0.005)+
  geom_line(data = data.frame(select(dt.GT,-method)),inherit.aes = FALSE,aes(x = X / N),stat="density",colour = "darkred",alpha = 0.3,size = 1.2)+
  geom_line(stat="density",alpha = 0.3,size = 1.2)+
  geom_vline(aes(xintercept = (0.5 + X.obs)/(0.5 + 0.5 + N),colour = method),alpha = 0.7,lty = "dashed")+
  theme_cowplot()


# math SE -----------------------------------------------------------------
library(magrittr)
library(extraDistr)
library(ggplot2)
library(cowplot)

set.seed(1)
N <- 1000L
n.rep <- 5000

alpha <- 200 + 0.5
beta <- N - alpha + 1

X0 <- rbinom(n.rep,N,alpha / (alpha + beta)) %>%
  data.frame(method = "X0",X = .)
X1 <- replicate(n.rep,rbeta(N,alpha,beta) %>% rbinom(N,1,.) %>% sum) %>%
  data.frame(method = "X1",X = .)
X2 <- replicate(n.rep,rbbinom(N,1,alpha,beta) %>% sum) %>%
  data.frame(method = "X2",X = .)
X3 <- rbbinom(n.rep,N,alpha,beta) %>%
  data.frame(method = "X3",X = .)
X4 <- replicate(n.rep,rbeta(1,alpha,beta) %>% rbinom(N,1,.) %>% sum) %>%
  data.frame(method = "X4",X = .)

rbind(X0,X1,X2,X3,X4) %>%
  ggplot(aes(X / N,fill = method,colour = method))+
  geom_density(alpha = .05)+
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,length = 5))+
  theme_minimal_grid()



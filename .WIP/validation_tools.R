# required packages -------------------------------------------------------
library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(extraDistr)
devtools::load_all()

# misc matrix tools ---------------------------------------------------------------------------
`upper.tri<-` <- function(x,value) {
  dx <- dim(x)
  if (length(dx) != 2L)
    stop("only matrix upper triangles can be replaced")
  len.i <- (prod(dx) - min(dx)) / 2
  len.v <- length(value)
  if (len.v != 1L && len.v != len.i)
    stop("replacement upper triangles has wrong length")
  x[upper.tri(x)] <- value
  x
}

`lower.tri<-` <- function(x,value) {
  dx <- dim(x)
  if (length(dx) != 2L)
    stop("only matrix lower triangles can be replaced")
  len.i <- (prod(dx) - min(dx)) / 2
  len.v <- length(value)
  if (len.v != 1L && len.v != len.i)
    stop("replacement lower triangles has wrong length")
  x[lower.tri(x)] <- value
  x
}

# network generating tools ------------------------------------------------

## generating P -----------------------------------------------------------

generate_P.seq <- function(n,node_names = as.character(1:n),
                           mode = c("directed","upper","lower","undirected")) {
  mode <- match.arg(mode)
  switch(mode,
         "directed" = {
           `subset_fun<-` <- SimuNet::`non.diagonal<-`
           l <- n^2 - n
         },
         "undirected" = ,
         "upper" = {
           `subset_fun<-` <- `upper.tri<-`
           l <- (n^2 - n) / 2
         },
         "lower" = {
           `subset_fun<-` <- `lower.tri<-`
           l <- (n^2 - n) / 2
         }
  )
  M <- matrix(0,n,n,dimnames = list(node_names,node_names))
  subset_fun(M) <-
    seq(0,1,length = l)
  M
}

generate_P.soc <- function(n,node_names = as.character(1:n)) {
  M <- matrix(0,n,n,dimnames = list(node_names,node_names))
  SimuNet::non.diagonal(M) <-
    seq(0,1,length = n^2 - n)
  M
}

generate_P.wrapper <- function(n,P.method = c("sequence","social"),node_names = as.character(1:n),...) {
  switch(P.method,
         "sequence" = generate_P.seq(n = n, node_names = node_names),
         "social" = generate_P.soc(n = n, node_names = node_names,...)
  )
}


## inferring P from observations ---------------------------------------
mat_rbeta <- function(A,N,alpha.prior = 0.5,beta.prior = 0.5) {
  P <- rbeta(n = length(A),A + alpha.prior,N - A + beta.prior)
  if (is.array(A)) {
    dim(P) <- dim(A)
    attr(P,"dimnames") <- attr(A,"dimnames")
    diag(P) <- 0
  }
  P
}

## draw scan(s) from P ----------------------------------------------------
mat_rbinom <- function(P,N) {
  Adj <- rbinom(length(P),N,P)
  if (is.array(P)) {
    dim(Adj) <- dim(P)
    attr(Adj,"dimnames") <- attr(P,"dimnames")
  }
  Adj
}

draw_scanList <- function(P = NULL,A = NULL,N,n.scans = N) {
  if (is.null(P) & !is.null(A)) {
    P <- mat_rbeta(A,N)
  }
  scan.list <- replicate(n = n.scans,simplify = FALSE,mat_rbinom(P,N = 1))
  class(scan.list) <- "sL"
  scan.list
}

print.sL <- function(x,...) {
  N <- length(x)
  if (N >= 10) {
    lapply(x[1:3],print)
    cat("\n... (",N - 3," more scans)\n\n",sep = "")
    lapply(x[(N - 2 + 1):N],print)
  } else {
    print.default(x)
  }
}


## draw adjacency from Beta-binomial  -------------------------------------
mat_rbbinom <- function(A,N,alpha.prior = 0.5,beta.prior = 0.5,N.new = N) {
  Adj <-
    extraDistr::rbbinom(
      n = length(A),
      size = N.new,
      alpha = A + alpha.prior,
      beta = N - A + beta.prior
    )
  if (is.array(A)) {
    dim(Adj) <- dim(A)
    attr(Adj,"dimnames") <- attr(A,"dimnames")
    diag(Adj) <- 0
  }
  Adj
}

## sum scans --------------------------------------------------------------
sum_scanList <- function(scan.list) {
  Reduce("+",scan.list)
}


## resample scanList -------------------------------------------------------
resample_scanList <- function(scan.list) {
  scan.list <- sample(scan.list,replace = TRUE)
  class(scan.list) <-  "sL"
  scan.list
}

bootstrap_scanList <- function(scan.list,n.boot = 100,simplify = FALSE) {
  replicate(n = n.boot,
            simplify = simplify,
            resample_scanList(scan.list) %>%
              sum_scanList()
  )
}


## wrapper for generating through each method -----------------------------
generate_infered_networks <- function(P,N,N.new = N,n.samp = 1000,seed = NULL) {
  if (!is.null(seed)) {set.seed(seed)}
  sL0 <- draw_scanList(P = P,n.scans = N)
  A0 <- sum_scan.list(sL0)

  list(
    A.GT = replicate(n = n.samp,simplify = FALSE,mat_rbinom(P = P,N = N)),
    A.bbinom = replicate(n = n.samp,simplify = FALSE,mat_rbbinom(A = A0,N = N)),
    A.SimuNet = replicate(n = n.samp,simplify = FALSE,mat_rbeta(A = A0,N = N) %>% mat_rbinom(P = .,N = N)),
    A.boot = bootstrap_scanList(scan.list = sL0,n.boot = n.samp)
  )
}

# data extraction tools ---------------------------------------------------

## extract xij ------------------------------------------------------------
extract_xij <- function(M,N = NULL,x.name = "x",...,
                        mode = c("directed","upper","lower","undirected"),keep.diag = FALSE) {
  n <- nrow(M)
  mode <- match.arg(mode)

  expand.grid(i = 1:n,j = 1:n) %>% {
    .$N <- N
    .[[x.name]] <- M[cbind(.$i,.$j)]
    .
  } %>%
    data.table() %>%
    cbind(...,.) %>%
    subset(switch(mode,"directed" = TRUE,"undirected" = ,"upper" = i <= j,"lower" = i >= j)) %>%
    subset(i != j | keep.diag)
}

create_Pij.dt <- function(P,mode = "directed") {
  P %>%
    extract_xij(x.name = "p",mode = mode,method = as.factor("original"))
}


extract_pij <- function(A.list,N,x.name = "p",keep.diag = FALSE,
                        mode = c("directed","upper","lower","undirected"),...) {
  lapply(A.list,function(A) A / N) %>%
    lapply(FUN = extract_xij,x.name = x.name,mode = mode,keep.diag = keep.diag,... = ...)

}

extract_pijs <- function(A,n,N,method = c("GT","bbinom","SimuNet","boot"),
                         x.name = "p",keep.diag = FALSE,
                         mode = c("directed","upper","lower","undirected"),...,n.samp) {
  mode <- match.arg(mode)
  n.edge <- switch(mode,
                   "directed" = n^2 - ifelse(keep.diag,0L,n),
                   "upper" = ,
                   "lower" = ,
                   "undirected" = (n^2 - ifelse(keep.diag,0L,n)) / 2
  )
  switch(method,
         "GT" = A$A.GT,
         "bbinom" = A$A.bbinom,
         "SimuNet" = A$A.SimuNet,
         "boot" = A$A.boot
  ) %>%
    extract_pij(N = N,x.name = x.name,mode = mode,keep.diag = keep.diag,...) %>%
    do.call(rbind,.) %>%
    cbind(
      n = n,
      N = N,
      method = factor(method),
      sample = factor(rep(1:n.samp,each = n.edge)),
      .
    )
}

extract_pijs.vec <- Vectorize(extract_pijs,
                              vectorize.args = "method",
                              SIMPLIFY = FALSE)


# calculate social network metrics ----------------------------------------

## network metric ---------------------------------------------------
extract_SNm <- function(Adj.list,SNm_fun,transform.to.igraph = FALSE,mode = NULL,weighted = TRUE,...,Xapply = lapply) {
  if (!is.list(Adj.list)) {Adj.list <- list(Adj.list)}
  if (transform.to.igraph) {
    Adj.list <- lapply(Adj.list,igraph::graph.adjacency,mode = mode,weighted = weighted)
  }
  Xapply(Adj.list,SNm_fun,...)
}


## wrapper to extract several metrics passed as a vector of functions ------


extract_SNm.vec <- Vectorize(FUN = extract_SNm,
                             vectorize.args = c("SNm_fun","transform.to.igraph"),
                             SIMPLIFY = FALSE)

### wrapper for node metrics ------------------------------------------------
custom_EV <- function(g) {
  igraph::eigen_centrality(g,directed = TRUE)$vector %>%
    data.table(EV = .)
}

custom_str <- function(g) {
  igraph::strength(g,mode = "all") %>%
    data.table(str = .)
}

custom_CC <- function(M) {
  c(DirectedClustering::ClustBCG(M,type = "directed")$totalCC) %>%
    data.table(CC = .)
}

custom_bet <- function(M) {
  sna::betweenness(M) %>%
    data.table(bet = .)
}

custom_fbet <- function(M) {
  sna::flowbet(M) %>%
    data.table(fbet = .)
}

#### wrappers for custom node metrics --------------------------------------------------
extract_node.metrics <- function(A,n,N,method = c("GT","bbinom","SimuNet","boot"),n.samp) {
  switch(method,
         "GT" = A$A.GT,
         "bbinom" = A$A.bbinom,
         "SimuNet" = A$A.SimuNet,
         "boot" = A$A.boot
  ) %>%
    lapply(function(A) A / N) %>%
    extract_SNm.vec(
      SNm_fun = c(
        custom_EV,
        custom_CC,
        custom_str,
        # custom_fbet,
        custom_bet
      ),
      transform.to.igraph = c(
        TRUE,
        FALSE,
        TRUE,
        # FALSE,
        FALSE
      ),
      mode = "directed",Xapply = rbind_lapply) %>%
    do.call(cbind,.) %>%
    cbind(
      n = n,
      N = N,
      method = factor(method),
      i = factor(1:n),
      sample = factor(rep(1:n.samp,each = n)),
      .
    )
}

extract_node.metrics.vec <- Vectorize(extract_node.metrics,
                                      vectorize.args = "method",
                                      SIMPLIFY = FALSE)

get_original.node.metrics <- function(P) {
  P %>%
  extract_SNm.vec(
    c(
      custom_EV,
      custom_CC,
      custom_str,
      # custom_fbet,
      custom_bet
    ),
    transform.to.igraph = c(
      TRUE,
      FALSE,
      TRUE,
      # FALSE,
      FALSE
    ),
    mode = "directed",Xapply = rbind_lapply) %>%
    do.call(cbind,.) %>%
    cbind(
      n = nrow(P),
      N = 1,
      method = factor("original"),
      i = factor(1:nrow(P)),
      sample = factor(1),
      .
    )
}

### wrappers for network metrics ------------
custom_GlobalCC <- function(M) {
  DirectedClustering::ClustBCG(M,type = "directed")$GlobaltotalCC %>%
    data.table(GCC = .)
}

custom_diam <- function(g) {
  igraph::diameter(g,directed = TRUE) %>%
    data.table(diam = .)
}

#### wrappers for custom global metrics --------------------------------------------------
extract_global.metrics <- function(A,n,N,method = c("GT","bbinom","SimuNet","boot"),n.samp) {
  switch(method,
         "GT" = A$A.GT,
         "bbinom" = A$A.bbinom,
         "SimuNet" = A$A.SimuNet,
         "boot" = A$A.boot
  ) %>%
    lapply(function(A) A / N) %>%
    extract_SNm.vec(SNm_fun = c(custom_GlobalCC,custom_diam),
                    transform.to.igraph = c(FALSE,TRUE),
                    mode = "directed",Xapply = rbind_lapply) %>%
    do.call(cbind,.) %>%
    cbind(
      n = n,
      N = N,
      method = factor(method),
      sample = factor(1:n.samp),
      .
    )
}

extract_global.metrics.vec <- Vectorize(extract_global.metrics,
                                        vectorize.args = "method",
                                        SIMPLIFY = FALSE)

get_original.global.metrics <- function(P) {
  P %>%
    extract_SNm.vec(c(custom_GlobalCC,custom_diam),
                    transform.to.igraph = c(FALSE,TRUE),
                    mode = "directed",Xapply = rbind_lapply) %>%
    do.call(cbind,.) %>%
    cbind(
      n = nrow(P),
      N = 1,
      method = factor("original"),
      sample = 1,
      .
    )
}

# summarize data -------------------------------------------------------------------------
## CI wrappers ----
CI_q <- function(x,probs = c(0.025,0.975)){
  quantile(x = x,probs = probs)
}

low_q <- function(x,probs = 0.025){
  quantile(x = x,probs = probs)
}

up_q <- function(x,probs = 0.975){
  quantile(x = x,probs = probs)
}

## data wrangling --------------
### general functions ----
calculate_CIs <- function(long.dt,by) {
  long.dt[,by = eval(substitute(by)),
               .(
                 value = median(value),
                 low = low_q(value),
                 up = up_q(value),
                 mean = mean(value),
                 sd = sd(value),
                 n.samp = .N
               )
  ]
}

compare_with_CI <- function(merged.CI.dt) {
  merged.CI.dt[,lower.CI := 0L,
  ][original < low,lower.CI := 1L,
  ][,higher.CI := 0L,
  ][up < original,higher.CI := 1L,
  ][,in.CI := 0L,
  ][lower.CI == 0L & higher.CI == 0L,in.CI := 1L,][]
}

proportion_in_CI <- function(compared.CI.dt,by) {
  compared.CI.dt[,by = eval(substitute(by)),
      .(
        in.CI = sum(.SD$in.CI),
        lower.CI = sum(.SD$lower.CI),
        higher.CI = sum(.SD$higher.CI),
        n.rep = .N
      )
    ][
      ,prop.in := in.CI / n.rep,
    ][
      ,prop.lower := lower.CI / n.rep,
    ][
      ,prop.higher := higher.CI / n.rep,
    ][]
}

### P_hat ----
calculate_CIs.P_hat <- function(P_hat.dt) {
  P_hat.dt[
    ,ij := factor(paste0(i,"-",j)),
  ][
    ,by = .(n,N,method,i,j,ij),
    .(
      p = median(p),
      low = low_q(p),
      up = up_q(p),
      mean = mean(p),
      sd = sd(p),
      n.samp = .N
    )
  ]
}

merge_with_original.P_hat <- function(P_hat.CI.dt,Pij.dt) {
  Pij.dt[,ij := factor(paste0(i,"-",j)),][]
  merge.data.table(
    x = P_hat.CI.dt,
    y = Pij.dt %>%
      select(ij,p) %>%
      setnames("p","original"),
    by = c("ij")
  )
}

proportion_in_CI.P_hat <- function(compared.CI.dt) {
  proportion_in_CI(compared.CI.dt,by = .(ij,n,N,method,i,j))
}

### nodes ----
melt_node.metrics <- function(node.metrics.dt) {
  node.metrics.dt %>%
    melt.data.table(
      id.vars = c("n","N","method","i","sample"),
      measure.vars = c(
        "EV",
        "CC",
        "str",
        # "fbet",
        "bet"
      ),
      variable.name = "metric",
      value.name = "value"
    )
}

calculate_CIs.nodes <- function(long.dt) {
  calculate_CIs(long.dt = long.dt,by = .(n,N,method,i,metric))
}

merge_with_original.nodes <- function(node.CI.dt,P.node.dt) {
  merge.data.table(
    x = node.CI.dt,
    y = P.node.dt %>%
      melt_node.metrics %>%
      select(i,metric,value) %>%
      setnames("value","original"),
    by = c("i","metric")
  )
}

proportion_in_CI.nodes <- function(compared.CI.dt) {
  proportion_in_CI(compared.CI.dt,by = .(n,N,i,method,metric))
}

### global ----
melt_global.metrics <- function(global.metrics.dt) {
  global.metrics.dt %>%
    melt.data.table(
      id.vars = c("n","N","method","sample"),
      measure.vars = c("GCC","diam"),
      variable.name = "metric",
      value.name = "value"
    )
}

calculate_CIs.global <- function(long.dt) {
  calculate_CIs(long.dt = long.dt,by = .(n,N,method,metric))
}

merge_with_original.global <- function(global.CI.dt,P.global.dt) {
  merge.data.table(
    x = global.CI.dt,
    y = P.global.dt %>%
      melt_global.metrics %>%
      select(metric,value) %>%
      setnames("value","original"),
    by = c("metric")
  )
}

proportion_in_CI.global <- function(compared.CI.dt) {
  proportion_in_CI(compared.CI.dt,by = .(n,N,method,metric))
}



# Repeat data inference for many networks ----------------------------------------------------
infer_one_network <- function(P,N,n.samp = n.samp,rep = NULL,Pij.dt,P.node.dt,P.global.dt) {
  n <- nrow(P)
  A0.list <- generate_infered_networks(P = P,N = N,n.samp = n.samp)
  P_hat.dt <-
    A0.list %>%
    extract_pijs.vec(method = c("GT","bbinom","SimuNet","boot"),n = n,N = N,n.samp = n.samp) %>%
    do.call(rbind,.) %>%
    calculate_CIs.P_hat %>%
    merge_with_original.P_hat(Pij.dt) %>%
    compare_with_CI %>%
    cbind(rep = rep,.)
  node_metric.dt <-
    A0.list %>%
    extract_node.metrics.vec(method = c("GT","bbinom","SimuNet","boot"),
                             n = n,N = N,n.samp = n.samp) %>%
    do.call(rbind,.) %>%
    melt_node.metrics() %>%
    calculate_CIs.nodes() %>%
    merge_with_original.nodes(P.node.dt) %>%
    compare_with_CI() %>%
    cbind(rep = rep,.)
  global_metric.dt <-
    A0.list %>%
    extract_global.metrics.vec(method = c("GT","bbinom","SimuNet","boot"),
                               n = n,N = N,n.samp = n.samp) %>%
    do.call(rbind,.)  %>%
    melt_global.metrics() %>%
    calculate_CIs.global() %>%
    merge_with_original.global(P.global.dt) %>%
    compare_with_CI() %>%
    cbind(rep = rep,.)
  list(P_hat.dt = P_hat.dt,node_metric.dt = node_metric.dt,global_metric.dt = global_metric.dt)
}

infer_multiple_networks <- function(P,n,N,n.rep,n.samp,cl = NULL,Pij.dt,P.node.dt,P.global.dt) {
  cat(paste0("\nn = ",n," - N = ",N," - n.rep = ",n.rep," - n.samp = ",n.samp,"\n"))
  pbapply::pblapply(1:n.rep,function(r) infer_one_network(P = P,N = N,n.samp = n.samp,rep = r,
                                                          Pij.dt = Pij.dt,P.node.dt = P.node.dt,
                                                          P.global.dt = P.global.dt),cl = cl) %>%
    {
      list(
        P_hat.dt = rbind_lapply(.,function(r) r$P_hat.dt),
        node_metric.dt = rbind_lapply(.,function(r) r$node_metric.dt),
        global_metric.dt = rbind_lapply(.,function(r) r$global_metric.dt)
      )
    }
}

infer_across_N <- Vectorize(infer_multiple_networks,vectorize.args = c("N"),SIMPLIFY = FALSE)

# Repeat data inference for different n ----------------------------------------------------
generate_P_and_infer <- function(n,N,n.rep,n.samp,cl = NULL) {
  P <- generate_P.seq(n = n,mode = "directed")
  Pij.dt <- create_Pij.dt(P,mode = "directed")
  P.node.dt <- get_original.node.metrics(P)
  P.global.dt <- get_original.global.metrics(P)
  infer_across_N(P = P,n = n,N = N,n.rep = n.rep,n.samp = n.samp,
                 cl = cl,
                 Pij.dt = Pij.dt,P.node.dt = P.node.dt,
                 P.global.dt = P.global.dt)
}

infer_across_n_N <- Vectorize(generate_P_and_infer,vectorize.args = c("n"),SIMPLIFY = FALSE)

combine_inferred_across_n_N <- function(dt.list) {
  list(
    P_hat.dt = rbind_lapply(
      dt.list,
      function(r.n) rbind_lapply(
        r.n,
        function(r) r$P_hat.dt
      )
    ),
    node_metric.dt = rbind_lapply(
      dt.list,
      function(r.n) rbind_lapply(
        r.n,
        function(r) r$node_metric.dt
      )
    ),
    global_metric.dt = rbind_lapply(
      dt.list,
      function(r.n) rbind_lapply(
        r.n,
        function(r) r$global_metric.dt
      )
    )
  )
}

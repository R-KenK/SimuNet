library(data.table)
library(igraph)
library(ggplot2)
library(dplyr)
library(pbapply)
devtools::load_all(".")

# custom functions -----
## Network generation algorithms ----
### Highest ranks choose a random number of partners among most valuables----
draw_npartner <- function(n, lambda = 2) {
  n.partner <- vector(mode = "integer",length = n)
  remain <- n
  for (i in 1:n) {
    n.partner[i] <- sample(0:remain,1,prob = dpois(0:remain,lambda = lambda))
    remain <- remain - n.partner[i]
  }
  n.partner
}

build_edgelist <- function(ind) {
  lapply(
    1:n,
    function(i) {
      n.partner <- ind[name == as.character(i)]$n.partner
      ind[name != as.character(i)][order(value,decreasing = TRUE)]$name |>
        head(n = n.partner) %>%
        cbind(i = as.character(i),j = .)
    }
  ) |>
    lapply(function(el) if (ncol(el) == 2) el else NULL) |>
    do.call(what = rbind)
}

highRankChooseMostValuable_generate_bin <- function(ind,samp.eff = 1) {
  replicate(n = samp.eff,
            {
              Adj <- matrix(0L,n,n,dimnames = list(ind$name,ind$name))
              ind$n.partner <- draw_npartner(n,2)[ind$rank]
              ind |>
                build_edgelist() %>%
                {Adj[.] <- 1L;Adj} %>%
                {ifelse(. + t(.) >= 1,1L,0L)} %>%
                {.[lower.tri(.)] <- 0L;.}
            }
  )
}

highRankChooseMostValuable_generate_asso <- function(samp.eff = 1,ind = NULL) {
  if (is.null(ind))
    ind <- data.table(name = as.character(1:n),value = runif(n),rank = sample(1:n,n))
  highRankChooseMostValuable_generate_bin(ind,samp.eff = samp.eff) |>
    rowSums(na.rm = TRUE,dims = 2L)
}

### Sum of ER graph (G(n,p = 0.5) ----
ER_generate_asso <- function(n,samp.eff,names = NULL,prob = 0.5) {
  M <- matrix(0L,n,n,dimnames = list(names,names))
  M[upper.tri(M)] <- rbinom(n = (n^2 - n) / 2,size = samp.eff,prob = prob)
  M
}

### Random graph (fixed p_ij from runif) ----
fixedRandom_generate_asso <- function(n,samp.eff,names = NULL,fixed.rand.unif = NULL) {
  if (is.null(fixed.rand.unif))
    fixed.rand.unif <- runif(n = (n^2 - n) / 2)
  M <- matrix(0L,n,n,dimnames = list(names,names))
  M[upper.tri(M)] <- rbinom(n = (n^2 - n) / 2,size = samp.eff,prob = fixed.rand.unif)
  M
}

### Random graph (varying p_ij from runif) ----
totalRandom_generate_asso <- function(n,samp.eff,names = NULL) {
  total.rand <- runif((n^2 - n) / 2)
  M <- matrix(0L,n,n,dimnames = list(names,names))
  M[upper.tri(M)] <- rbinom(n = (n^2 - n) / 2,size = samp.eff,prob = total.rand)
  M
}

## Network data extraction ----
compute.EV <- function(M) {
  M |>
    igraph::graph.adjacency(mode = "upper",weighted = TRUE)  %>%
    {igraph::eigen_centrality(.)$vector}
}

## Network comparison ----
measure_distances <- function(edge.dt,x,y) {
  edge.comp.dt <-
    edge.dt[
      ,
      list(
        x = list(.SD[type == x]$weight),
        y = list(.SD[type == y]$weight)
      ),
      by = .(i,j,n,samp.eff,n.rep),.SDcols = "weight"
    ]
  edge.comp.dt[,possible := lapply(samp.eff,\(r) 0:r),]

  edge.comp.dt[,x.table := lapply(1:nrow(edge.comp.dt),\(r) factor(x[[r]],levels = possible[[r]]) |> table()),]
  edge.comp.dt[,y.table := lapply(1:nrow(edge.comp.dt),\(r) factor(y[[r]],levels = possible[[r]]) |> table()),]
  edge.comp.dt[,x.prob := lapply(1:nrow(edge.comp.dt),\(r) x.table[[r]] / sum(x.table[[r]])),]
  edge.comp.dt[,y.prob := lapply(1:nrow(edge.comp.dt),\(r) y.table[[r]] / sum(y.table[[r]])),]
  edge.comp.dt[,x.mat := lapply(seq_along(x.table),\(r) cbind(possible[[r]],x.table[[r]])),]
  edge.comp.dt[,y.mat := lapply(seq_along(y.table),\(r) cbind(possible[[r]],y.table[[r]])),]
  edge.comp.dt[,x.mean := sapply(x,mean),]
  edge.comp.dt[,y.mean := sapply(y,mean),]
  edge.comp.dt[,mean.diff := sapply(1:nrow(edge.comp.dt),\(r) abs(x.mean[[r]] - y.mean[[r]])),]
  edge.comp.dt[,KS.stat := sapply(1:nrow(edge.comp.dt),\(r) ks.test(x[[r]],y[[r]])$statistic),]
  edge.comp.dt[,KS.p := sapply(1:nrow(edge.comp.dt),\(r) ks.test(x[[r]],y[[r]])$p.value),]
  edge.comp.dt[,KL := sapply(1:nrow(edge.comp.dt),\(r) philentropy::KL(rbind(x.prob[[r]],y.prob[[r]]))),]
  edge.comp.dt[,JS := sapply(1:nrow(edge.comp.dt),\(r) philentropy::JSD(rbind(x.prob[[r]],y.prob[[r]]))),]
  edge.comp.dt[,EMD.e := sapply(seq_along(x.table),\(r) emdist::emd(x.mat[[r]],y.mat[[r]],dist = "euclidean")),]
  edge.comp.dt[,EMD.m := sapply(seq_along(x.table),\(r) emdist::emd(x.mat[[r]],y.mat[[r]],dist = "manhattan")),]

  edge.comp.dt |>
    select(c("i","j","n","samp.eff","n.rep","mean.diff","KS.stat","KS.p","KL","JS","EMD.e","EMD.m")) %>%
    cbind(x = x,y = y,.)
}

## Encapsulating case study ----
run_experiment <- function(netgen_fun,netgen_args,netgen_other,
                           n,samp.eff,names = NULL,n.rep,
                           cl,ind = NULL,ind.other = NULL,fixed.rand.prob = NULL) {
  if (is.null(fixed.rand.prob)) {
    fixed.rand.prob <- runif(n = (n^2 - n) / 2)
  }
  if (is.null(ind) | is.null(ind.other)) {
    if (is.null(ind)) {
      cat("Generating new individuals...\n")
      ind <- data.table(
        name = as.character(1:n),
        value = runif(n),
        rank = sample(1:n,n)
      )
    }
    if (is.null(ind.other)) {
      cat("Generating new individuals (other group)...\n")
      ind.other <- data.table(
        name = as.character(1:n),
        value = runif(n),
        rank = sample(1:n,n)
      )
    }
    cat("Exporting new individuals to clusters...\n")
  } else {
    cat("Using provided individuals...\n")
  }

  snow::clusterExport(
    cl,
    list("n","samp.eff","ind","ind.other","fixed.rand.prob"),
    envir = environment()
  )

  cat("Running the simulations...\n")
  cat(" + Generating many series of associations for the real network:\n")
  real.dist <- pbapply::pbreplicate(
    n = n.rep,cl = cl,
    expr = do.call(netgen_fun,netgen_args)
  )
  cat(" + Generating many series of associations for another network:\n")
  other.dist <-  pbapply::pbreplicate(
    n = n.rep,cl = cl,
    expr = do.call(netgen_fun,netgen_other)
  )
  cat(" + Generating many weighted random graphs (equivalent to the sum of G(n,p = 0.5)):\n")
  ER.dist <- pbapply::pbreplicate(
    n = n.rep,cl = cl,
    expr = ER_generate_asso(n = n,samp.eff = samp.eff,names = names)
  )
  cat(" + Generating many weighted random graphs (with random but fixed p_ij):\n")
  fixed.rand.dist <- pbapply::pbreplicate(
    n = n.rep,cl = cl,
    expr = fixedRandom_generate_asso(n = n,samp.eff = samp.eff,
                                     names = names,fixed.rand.unif = fixed.rand.prob)
  )
  cat(" + Generating many weighted random graphs (with random p_ij varying at each replication):\n")
  total.rand.dist <- pbapply::pbreplicate(
    n = n.rep,cl = cl,
    expr = totalRandom_generate_asso(n = n,samp.eff = samp.eff,names = names)
  )
  cat(" + Generating networks with SimuNet from the weighted adjacency matrix observed:\n")
  SimuNet.dist <- pbapply::pbreplicate(
    n = n.rep,cl = cl,
    expr = real.dist[,,1] |>
      simunet(samp.effort = samp.eff,mode = "upper",n.scans = samp.eff,
              alpha.prior = 1,beta.prior = 1) |>
      sum_scans()
  )

  cat("Aggregating data...\n")
  edge.dt <-
    expand.grid(i = 1:n,j = 1:n,rep = 1:n.rep) |>
    subset(i < j) %>%
    {
      rbind(
        cbind(.,n,samp.eff = samp.eff,n.rep = n.rep,type = "real",
              weight = real.dist[as.matrix(.)]),
        cbind(.,n,samp.eff = samp.eff,n.rep = n.rep,type = "other",
              weight = other.dist[as.matrix(.)]),
        cbind(.,n,samp.eff = samp.eff,n.rep = n.rep,type = "ER",
              weight = ER.dist[as.matrix(.)]),
        cbind(.,n,samp.eff = samp.eff,n.rep = n.rep,type = "fixed.rand",
              weight = fixed.rand.dist[as.matrix(.)]),
        cbind(.,n,samp.eff = samp.eff,n.rep = n.rep,type = "total.rand",
              weight = total.rand.dist[as.matrix(.)]),
        cbind(.,n,samp.eff = samp.eff,n.rep = n.rep,type = "SimuNet",
              weight = SimuNet.dist[as.matrix(.)])
      )
    } |>
    data.table()

  cat("Measuring distances...\n")
  edge.comp.dt <-
    rbind(
      suppressWarnings(suppressMessages(measure_distances(edge.dt,"SimuNet","real"))),
      suppressWarnings(suppressMessages(measure_distances(edge.dt,"SimuNet","other"))),
      suppressWarnings(suppressMessages(measure_distances(edge.dt,"SimuNet","ER"))),
      suppressWarnings(suppressMessages(measure_distances(edge.dt,"SimuNet","fixed.rand"))),
      suppressWarnings(suppressMessages(measure_distances(edge.dt,"SimuNet","total.rand")))
    )
  edge.comp.dt$y <- edge.comp.dt$y |> factor(levels = c("real","other","ER","fixed.rand","total.rand"))
  edge.dt$type <- edge.dt$type |> factor(levels = c("SimuNet","real","other","ER","fixed.rand","total.rand"))
  list(real.dist = real.dist,
       other.dist = other.dist,
       ER.dist = ER.dist,
       fixed.rand.dist = fixed.rand.dist,
       total.rand.dist = total.rand.dist,
       SimuNet.dist = SimuNet.dist,
       edge.dt = edge.dt,
       edge.comp.dt = edge.comp.dt)
}

## Aggregating experiment results ----
rbind_2results <- function(res.1,res.2) {
  res <- lapply(
    seq_along(res.1),
    function(l) {
      if (is.data.frame(res.1[[l]]))
        bind_fun <- rbind
      else
        bind_fun <- function(...) Reduce(function(x,y) array(c(x,y),dim = c(dim(x)[1:2],dim(x)[3] + dim(y)[3])),list(...))
      bind_fun(
        res.1[[l]],
        res.2[[l]]
      )
    }
  )
  names(res) <- names(res.1)
  res
}

rbind_results <- function(...) {
  Reduce(rbind_2results,list(...))

}

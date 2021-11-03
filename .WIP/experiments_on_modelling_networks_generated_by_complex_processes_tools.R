library(data.table)
library(igraph)
library(ggplot2)
library(dplyr)
library(pbapply)
library(ggridges)
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

build_edgelist <- function(group) {
  lapply(
    1:n,
    function(i) {
      n.partner <- group[name == as.character(i)]$n.partner
      group[name != as.character(i)][order(value,decreasing = TRUE)]$name |>
        head(n = n.partner) %>%
        cbind(i = as.character(i),j = .)
    }
  ) |>
    lapply(function(el) if (ncol(el) == 2) el else NULL) |>
    do.call(what = rbind)
}

highRankChooseMostValuable_generate_bin <- function(group,samp.eff = 1) {
  replicate(n = samp.eff,
            {
              Adj <- matrix(0L,n,n,dimnames = list(group$name,group$name))
              group$n.partner <- draw_npartner(n,2)[group$rank]
              group |>
                build_edgelist() %>%
                {Adj[.] <- 1L;Adj} %>%
                {ifelse(. + t(.) >= 1,1L,0L)} %>%
                {.[lower.tri(.)] <- 0L;.}
            }
  )
}

highRankChooseMostValuable_generate_asso <- function(samp.eff = 1,group = NULL) {
  if (is.null(group))
    group <- data.table(name = as.character(1:n),value = runif(n),rank = sample(1:n,n))
  highRankChooseMostValuable_generate_bin(group,samp.eff = samp.eff) |>
    rowSums(na.rm = TRUE,dims = 2L)
}

### Most gregarious associate more likely, within cliques ----
greg_product <- function(group) {
  M <- group$gregariousness %o% group$gregariousness
  rownames(M) <- colnames(M) <- group$name
  M[lower.tri(M,diag = TRUE)] <- 0
  M
}

is_sameClique <- function(group) {
  expand.grid(i = 1:nrow(group),j = 1:nrow(group)) |>
    subset(i < j) %>%
    {.$same.clique <- group$clique[.$i] == group$clique[.$j];.} %>% {
      M <- matrix(NA,nrow(group),nrow(group),
                  dimnames = list(as.character(1:nrow(group)),as.character(1:nrow(group)))
      )
      M[cbind(.$i,.$j)] <- .$same
      M
    }
}

draw_matbinom <- function(M,samp.eff) {
  if (all(M == 0)) return(M)
  M[upper.tri(M)] <-
    rbinom(length(M[upper.tri(M)]),size = samp.eff,prob = M[upper.tri(M)])
  M
}

gregWithinClique_generate_asso <- function(group,samp.eff,...) {
  same.clique <- is_sameClique(group = group)
  within <-
    ifelse(same.clique & !is.na(same.clique),greg_product(group),0)
  out <-
    ifelse(!same.clique & !is.na(same.clique),greg_product(group) / 10,0)
  {within + out} |>
    draw_matbinom(samp.eff = samp.eff)
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

## Parameter list generation ----
generate_paramList <- function(param.n,param.samp.eff,param.netgen,
                               n.rep,n.group,n.each) {
  # data.table expandgrid equivalent (handles lists too, not only vectors)
  param.list <-
    CJ(
      n = param.n,
      samp.eff = param.samp.eff,
      netgen = param.netgen,
      n.rep = n.rep,
      group.number = 1:n.group,
      group.rep = 1:n.each,
      sorted = FALSE
    )

  # "extracts" from the list of names & functions into their own columns
  param.list[,netgen_name := sapply(.SD$netgen,"[",j = "netgen_name") |> unlist(),]
  param.list[,netgen_fun := sapply(.SD$netgen,\(dt) dt[j = "netgen_fun"][[1]]),]
  param.list[,netgen := NULL,]
  param.list

  param.list <-
    # generate one group per group.number
    param.list[,by = .(n,samp.eff,netgen_name,n.rep,group.number),
               .(netgen_args = list(
                 group = data.table(
                   name = as.character(1:n),
                   value = runif(n),
                   rank = sample(1:n,n),
                   gregariousness = rbeta(n,0.5,0.5),
                   clique = sample(LETTERS[1:(rpois(1,1) + 1)],n,replace = TRUE)
                 )
               )
               )
    ] |>
    merge.data.table(x = param.list,by = c("n","samp.eff","netgen_name","n.rep","group.number")) %>%
    # generate one other group per group.number & group.rep
    {
      .[,by = .(n,samp.eff,netgen_name,n.rep,group.number,group.rep),
        .(
          netgen_other = list(
            group = data.table(
              name = as.character(1:n),
              value = runif(n),
              rank = sample(1:n,n),
              gregariousness = runif(n),
              clique = sample(LETTERS[1:(rpois(1,1) + 1)],n,replace = TRUE)
            )
          ),
          fixed.rand.prob = list(runif(n = (n^2 - n) / 2))
        )
      ] |>
        merge.data.table(x = .,by = c("n","samp.eff","netgen_name","n.rep","group.number","group.rep"))
    } %>%
    # aggregates into a list the arguments to be passed to the netgen_fun
    .[,
      `:=`(netgen_args = mapply(list,samp.eff = samp.eff,group = netgen_args,
                                fixed.rand.prob = fixed.rand.prob,SIMPLIFY = FALSE),
           netgen_other = mapply(list,samp.eff = samp.eff,group = netgen_other,SIMPLIFY = FALSE)
      )
    ] %>%
    .[,fixed.rand.prob := NULL]
  param.list
}

## Network data extraction ----
compute.EV <- function(M) {
  M |>
    igraph::graph.adjacency(mode = "upper",weighted = TRUE)  %>%
    {igraph::eigen_centrality(.)$vector}
}

## Network comparison ----
prepare_for_distances <- function(param.list) {
  message("Calculating weight probabilities...")
  param.list[,possible    := lapply(samp.eff,\(r) 0:r)]
  param.list[,edge.weight := lapply(edge.dt,\(dt) dt[,c("type","i","j","weight")])]
  param.list[,edge.weight := lapply(edge.weight,\(dt) dt[,.(weight = list(.SD$weight)),by = .(type,i,j)])]
  param.list[,edge.weight := mapply(
    \(ew.dt,pos) ew.dt[
      ,
      weight.prob := lapply(weight,factor,levels = pos) |>
        lapply(table) |>
        lapply(\(w) w / sum(w))
    ][
      ,
      weight.mat := lapply(weight.prob,\(prob) cbind(pos,prob))
    ],
    ew.dt = edge.weight,pos = possible,SIMPLIFY = FALSE)]
  message("Formatting data...")
  param.list[,edge.weight := lapply(edge.weight,\(ew.dt) ew.dt[,weight.mean := sapply(weight,mean)])]
  param.list[
    ,
    edge.weight :=
      lapply(
        edge.weight,
        \(ew.dt) lapply(unique(ew.dt$type),\(ty) ew.dt[type == ty]) %>%
          {names(.) <- unique(ew.dt$type);.}
      )
  ]
  param.list[]
}

measure_distances <- function(param.list,x,y) {
  message("Measuring distances ",x,"-",y,"...")
  edge.weight.distance <-
    lapply(
      param.list$edge.weight,
      \(ew.list) cbind(
        reference = ew.list[[x]]$type,
        ew.list[[y]][,c("type","i","j")],
        data.table(
          mean.diff = abs(ew.list[[y]]$weight.mean - ew.list[[x]]$weight.mean),
          KS.stat   = mapply(\(x,y) suppressWarnings(ks.test(x,y,exact = FALSE)$statistic),
                             x = ew.list[[x]]$weight,y = ew.list[[y]]$weight),
          KS.p      = mapply(\(x,y) suppressWarnings(ks.test(x,y,exact = FALSE)$p.value),
                             x = ew.list[[x]]$weight,y = ew.list[[y]]$weight),
          KL        = mapply(\(x,y) suppressMessages(philentropy::KL(rbind(x,y))),
                             x = ew.list[[x]]$weight.prob,y = ew.list[[y]]$weight.prob),
          JS        = mapply(\(x,y) suppressMessages(philentropy::JSD(rbind(x,y))),
                             x = ew.list[[x]]$weight.prob,y = ew.list[[y]]$weight.prob),
          EMD.e     = mapply(\(x,y) emdist::emd(x,y,dist = "euclidean"),
                             x = ew.list[[x]]$weight.mat,y = ew.list[[y]]$weight.mat),
          EMD.m     = mapply(\(x,y) emdist::emd(x,y,dist = "manhattan"),
                             x = ew.list[[x]]$weight.mat,y = ew.list[[y]]$weight.mat)
        )
      )
    )
  mapply(cbind,n = param.list$n,samp.eff = param.list$samp.eff,n.rep = param.list$n.rep,
         group.number = param.list$group.number,group.rep = param.list$group.rep,
         edge.weight.distance,SIMPLIFY = FALSE
  )
}

measure_all_distances <- function(param.list,x = "real",y = c("real.bis","other","ER",
                                                              "fixed.rand","total.rand","SimuNet"),
                                  cl) {
  distances <- data.table(x = x,y = y)
  param.list$edge.distances.dt <-
    pbapply::pblapply(
      1:nrow(distances),
      \(r) {
        x <- distances$x[r];y <- distances$y[r]
        measure_distances(param.list,x,y)
      },cl = cl
    ) |>
    split_seq() |>
    lapply(do.call,what = rbind)
  param.list[,":="(possible = NULL,edge.weight = NULL)]
  param.list[]
}


## Encapsulating case study ----
# run_experiment <- function(netgen_fun,netgen_args,netgen_other,
#                            n,samp.eff,names = NULL,n.rep,
#                            cl,group = NULL,group.other = NULL,fixed.rand.prob = NULL) {
#   if (is.null(fixed.rand.prob)) {
#     fixed.rand.prob <- runif(n = (n^2 - n) / 2)
#   }
#   if (is.null(group) | is.null(group.other)) {
#     if (is.null(group)) {
#       message("Generating new individuals...")
#       group <- data.table(
#         name = as.character(1:n),
#         value = runif(n),
#         rank = sample(1:n,n)
#       )
#     }
#     if (is.null(group.other)) {
#       message("Generating new individuals (other group)...")
#       group.other <- data.table(
#         name = as.character(1:n),
#         value = runif(n),
#         rank = sample(1:n,n)
#       )
#     }
#     message("Exporting new individuals to clusters...")
#   } else {
#     message("Using provided individuals...")
#   }
#
#   snow::clusterExport(
#     cl,
#     list("n","samp.eff","group","group.other","fixed.rand.prob"),
#     envir = environment()
#   )
#
#   message("Running the simulations...")
#   message(" + Generating many series of associations for the real network:")
#   real.dist.double <- pbapply::pbreplicate(
#     n = n.rep * 2,cl = cl,
#     expr = do.call(netgen_fun,netgen_args)
#   )
#   real.dist <- real.dist.double[,,1:n.rep]
#   real.dist.bis <- real.dist.double[,,(n.rep + 1):dim(real.dist.double)[3]]
#   message(" + Generating many series of associations for another network:")
#   other.dist <-  pbapply::pbreplicate(
#     n = n.rep,cl = cl,
#     expr = do.call(netgen_fun,netgen_other)
#   )
#   message(" + Generating many weighted random graphs (equivalent to the sum of G(n,p = 0.5)):")
#   ER.dist <- pbapply::pbreplicate(
#     n = n.rep,cl = cl,
#     expr = ER_generate_asso(n = n,samp.eff = samp.eff,names = names)
#   )
#   message(" + Generating many weighted random graphs (with random but fixed p_ij):")
#   fixed.rand.dist <- pbapply::pbreplicate(
#     n = n.rep,cl = cl,
#     expr = fixedRandom_generate_asso(n = n,samp.eff = samp.eff,
#                                      names = names,fixed.rand.unif = fixed.rand.prob)
#   )
#   message(" + Generating many weighted random graphs (with random p_ij varying at each replication):")
#   total.rand.dist <- pbapply::pbreplicate(
#     n = n.rep,cl = cl,
#     expr = totalRandom_generate_asso(n = n,samp.eff = samp.eff,names = names)
#   )
#   message(" + Generating networks with SimuNet from the weighted adjacency matrix observed:")
#   SimuNet.dist <- pbapply::pbreplicate(
#     n = n.rep,cl = cl,
#     expr = real.dist[,,1] |>
#       simunet(samp.effort = samp.eff,mode = "upper",n.scans = samp.eff,
#               alpha.prior = 1,beta.prior = 1) |>
#       sum_scans()
#   )
#
#   message("Aggregating data...")
#   edge.dt <-
#     expand.grid(i = 1:n,j = 1:n,rep = 1:n.rep) |>
#     subset(i < j) %>%
#     {
#       rbind(
#         cbind(.,n,samp.eff = samp.eff,n.rep = n.rep,type = "real",
#               weight = real.dist[as.matrix(.)]),
#         cbind(.,n,samp.eff = samp.eff,n.rep = n.rep,type = "real.bis",
#               weight = real.dist.bis[as.matrix(.)]),
#         cbind(.,n,samp.eff = samp.eff,n.rep = n.rep,type = "other",
#               weight = other.dist[as.matrix(.)]),
#         cbind(.,n,samp.eff = samp.eff,n.rep = n.rep,type = "ER",
#               weight = ER.dist[as.matrix(.)]),
#         cbind(.,n,samp.eff = samp.eff,n.rep = n.rep,type = "fixed.rand",
#               weight = fixed.rand.dist[as.matrix(.)]),
#         cbind(.,n,samp.eff = samp.eff,n.rep = n.rep,type = "total.rand",
#               weight = total.rand.dist[as.matrix(.)]),
#         cbind(.,n,samp.eff = samp.eff,n.rep = n.rep,type = "SimuNet",
#               weight = SimuNet.dist[as.matrix(.)])
#       )
#     } |>
#     data.table()
#
#   message("Measuring distances...")
#   edge.distance.dt <-
#     rbind(
#       suppressWarnings(suppressMessages(measure_distances(edge.dt,"real","real.bis"))),
#       suppressWarnings(suppressMessages(measure_distances(edge.dt,"SimuNet","real"))),
#       suppressWarnings(suppressMessages(measure_distances(edge.dt,"SimuNet","other"))),
#       suppressWarnings(suppressMessages(measure_distances(edge.dt,"SimuNet","ER"))),
#       suppressWarnings(suppressMessages(measure_distances(edge.dt,"SimuNet","fixed.rand"))),
#       suppressWarnings(suppressMessages(measure_distances(edge.dt,"SimuNet","total.rand")))
#     )
#   edge.distance.dt$y <- edge.distance.dt$y |> factor(levels = c("real.bis","real","other","ER","fixed.rand","total.rand"))
#   edge.dt$type <- edge.dt$type |> factor(levels = c("SimuNet","real","real.bis","other","ER","fixed.rand","total.rand"))
#   list(real.dist = real.dist,
#        real.dist.bis = real.dist.bis,
#        other.dist = other.dist,
#        ER.dist = ER.dist,
#        fixed.rand.dist = fixed.rand.dist,
#        total.rand.dist = total.rand.dist,
#        SimuNet.dist = SimuNet.dist,
#        edge.dt = edge.dt,
#        edge.distance.dt = edge.distance.dt)
# }
## Encapsulating simulations ----
# run_simulations <- function(param.list,verbose = TRUE) {
#   lapply(
#     1:nrow(param.list),
#     function(r) {
#       n <- param.list$n[r]; samp.eff <- param.list$samp.eff[r]; n.rep <- param.list$n.rep[r]
#       group <- param.list$group[[r]];netget_fun <- param.list$netgen_fun[[r]]
#       group.number <- param.list$group.number[[r]]; group.rep <- param.list$group.rep[r]
#       if (verbose) {
#         message(paste0("\nn = ",n," - samp.eff = ",samp.eff," - group.rep = ",group.rep,
#                        " (params ",r,"/",nrow(param.list),")"))
#         message("________________________________")
#         results <- run_experiment(
#           netgen_fun = highRankChooseMostValuable_generate_asso,
#           netgen_args = list(samp.eff = samp.eff,group = group),
#           netgen_other = list(samp.eff = samp.eff,group = group.other),
#           n = n,samp.eff = samp.eff,names = group$name,n.rep = n.rep,cl = cl
#         )
#       } else {
#         message(paste0("n = ",n," - samp.eff = ",samp.eff," - group.rep = ",group.rep,
#                        " (params ",r,"/",nrow(param.list),")    - ",Sys.time()))
#         results <- purrr::quietly(run_experiment)(
#           netgen_fun = highRankChooseMostValuable_generate_asso,
#           netgen_args = list(samp.eff = samp.eff,group = group),
#           netgen_other = list(samp.eff = samp.eff,group = group.other),
#           n = n,samp.eff = samp.eff,names = group$name,n.rep = n.rep,cl = cl
#         )$result
#       }
#       results$edge.dt <- cbind(results$edge.dt,group.number = group.number,group.rep = group.rep)
#       results$edge.distance.dt <- cbind(results$edge.distance.dt,group.number = group.number,group.rep = group.rep)
#       if (verbose) message("________________________________\n\n")
#       results
#     }
#   )
# }
split_into_chunks <- function(x,n) {
  chunks <- split(x, sort(x%%n))
  names(chunks) <- NULL
  chunks
}

split_seq <- function(l) {
  n <- length(l[[1]])
  l <- split(unlist(l, use.names = FALSE,recursive = FALSE), seq_len(n))
  names(l) <- NULL
  l
}


prepare_clusters <- function(n.cores = 7L,param.list) {
  message("Creating parallel clusters...")
  cl <- snow::makeCluster(7)
  snow::clusterCall(
    cl,function() {
      source(".WIP/experiments_on_modelling_networks_generated_by_complex_processes_tools.R")
    }
  )
  snow::clusterExport(cl,list("param.list"),envir = environment())
  cl
}

run_simulations <- function(param.list,n.steps = 5L) {
  message("Running the simulations...")

  if (is.null(cl))
    n.chunks <- 1L
  else
    n.chunks <- length(cl) * n.steps

  param.list$netgen_output <-
    pbapply::pblapply(
      split_into_chunks(1:nrow(param.list),n.chunks),
      function(chunk) {
        lapply(
          chunk,
          \(r) {
            n <- param.list$n[r];
            samp.eff <- param.list$samp.eff[r];
            n.rep <- param.list$n.rep[r];
            names <- param.list$netgen_args[[r]]$group$name;
            fixed.rand.prob <- param.list$netgen_args[[r]]$fixed.rand.prob;
            f <- param.list$netgen_fun[[r]];
            args <- param.list$netgen_args[[r]];
            other <- param.list$netgen_other[[r]];

            real <- replicate(
              n = n.rep,
              expr = do.call(what = f,args = args) |>
                {\(x) {mode(x) <- "integer";x}}()
            )
            real.bis <- replicate(
              n = n.rep,
              expr = do.call(what = f,args = args) |>
                {\(x) {mode(x) <- "integer";x}}()
            )
            other <- replicate(
              n = n.rep,
              expr = do.call(what = f,args = other) |>
                {\(x) {mode(x) <- "integer";x}}()
            )
            ER <- replicate(
              n = n.rep,
              expr = ER_generate_asso(n = n,samp.eff = samp.eff,names = names)
            )
            fixed.rand <- replicate(
              n = n.rep,
              expr = fixedRandom_generate_asso(n = n,samp.eff = samp.eff,
                                               names = names,fixed.rand.unif = fixed.rand.prob)
            )
            total.rand <- replicate(
              n = n.rep,
              expr = totalRandom_generate_asso(n = n,samp.eff = samp.eff,
                                               names = names)
            )
            SimuNet <- replicate(
              n = n.rep,
              expr = real[,,1] |>
                simunet(samp.effort = samp.eff,mode = "upper",n.scans = samp.eff,
                        alpha.prior = 1,beta.prior = 1) |>
                sum_scans()
            )
            data.table(
              type = c("real",
                       "real.bis",
                       "other",
                       "ER",
                       "fixed.rand",
                       "total.rand",
                       "SimuNet"),
              dist = lapply(
                list(real,
                     real.bis,
                     other,
                     ER,
                     fixed.rand,
                     total.rand,
                     SimuNet),
                list
              )
            )
          }
        )
      },cl = cl
    ) |>
    unlist(recursive = FALSE)
  param.list
}

aggregate_edgeDT <- function(param.list,cl = NULL,n.steps = 2L) {
  if (is.null(cl))
    n.chunks <- 1L
  else
    n.chunks <- length(cl) * n.steps

  message("Aggregating edge weights...")

  param.list$edge.dt <-
    pbapply::pblapply(
      split_into_chunks(1:nrow(param.list),n.chunks),
      function(chunk) {
        lapply(
          chunk,
          \(r) {
            n <- param.list$n[r]
            samp.eff <- param.list$samp.eff[r]
            n.rep <- param.list$n.rep[r]
            group.number <- param.list$group.number[r]
            group.rep <- param.list$group.rep[r]
            real <- param.list$netgen_output[[r]][type == "real"]$dist[[1]][[1]]
            real.bis <- param.list$netgen_output[[r]][type == "real.bis"]$dist[[1]][[1]]
            other <- param.list$netgen_output[[r]][type == "other"]$dist[[1]][[1]]
            ER <- param.list$netgen_output[[r]][type == "ER"]$dist[[1]][[1]]
            fixed.rand <- param.list$netgen_output[[r]][type == "fixed.rand"]$dist[[1]][[1]]
            total.rand <- param.list$netgen_output[[r]][type == "total.rand"]$dist[[1]][[1]]
            SimuNet <- param.list$netgen_output[[r]][type == "SimuNet"]$dist[[1]][[1]]
            edge.dt <-
              expand.grid(i = 1:n,j = 1:n,rep = 1:n.rep) |>
              subset(i < j) %>%
              {
                rbind(
                  cbind(.,n = n,samp.eff = samp.eff,n.rep = n.rep,type = "real",
                        weight = real[as.matrix(.)]),
                  cbind(.,n = n,samp.eff = samp.eff,n.rep = n.rep,type = "real.bis",
                        weight = real.bis[as.matrix(.)]),
                  cbind(.,n = n,samp.eff = samp.eff,n.rep = n.rep,type = "other",
                        weight = other[as.matrix(.)]),
                  cbind(.,n = n,samp.eff = samp.eff,n.rep = n.rep,type = "ER",
                        weight = ER[as.matrix(.)]),
                  cbind(.,n = n,samp.eff = samp.eff,n.rep = n.rep,type = "fixed.rand",
                        weight = fixed.rand[as.matrix(.)]),
                  cbind(.,n = n,samp.eff = samp.eff,n.rep = n.rep,type = "total.rand",
                        weight = total.rand[as.matrix(.)]),
                  cbind(.,n = n,samp.eff = samp.eff,n.rep = n.rep,type = "SimuNet",
                        weight = SimuNet[as.matrix(.)])
                )
              } |>
              data.table()
            edge.dt$type <-
              edge.dt$type |>
              factor(levels = c("real","real.bis","SimuNet","other","ER","fixed.rand","total.rand"))
            cbind(edge.dt,group.number = group.number,group.rep = group.rep)
          }
        )
      },cl = cl
    ) |>
    unlist(recursive = FALSE)
  param.list
}

calculate_edgeDistanceDT <- function(param.list,cl) {
  param.list <- prepare_for_distances(param.list)
  measure_all_distances(param.list,cl = cl)
}

## Aggregating experiment results ----
# aggregate_results <- function(param.list,results.list,save.dists = FALSE) {
#   lapply(
#     1:nrow(param.list),
#     function(r) {
#       cbind(
#         param.list[r,],
#         data.table(
#           dists = list(data.table(
#             real.dist = list(results.list[[r]]$real.dist),
#             real.dist.bis = list(results.list[[r]]$real.dist.bis),
#             other.dist = list(results.list[[r]]$other.dist),
#             ER.dist = list(results.list[[r]]$ER.dist),
#             fixed.rand.dist = list(results.list[[r]]$fixed.rand.dist),
#             total.rand.dist = list(results.list[[r]]$total.rand.dist),
#             SimuNet.dist = list(results.list[[r]]$SimuNet.dist)
#           )),
#           edge.dt = list(results.list[[r]]$edge.dt),
#           edge.distance.dt = list(results.list[[r]]$edge.distance.dt)
#         )
#       )
#     }
#   ) |>
#     rbindlist()
# }

# rbind_2results <- function(res.1,res.2) {
#   res <- lapply(
#     seq_along(res.1),
#     function(l) {
#       if (is.data.frame(res.1[[l]]))
#         bind_fun <- rbind
#       else
#         bind_fun <- function(...) Reduce(function(x,y) array(c(x,y),dim = c(dim(x)[1:2],dim(x)[3] + dim(y)[3])),list(...))
#       bind_fun(
#         res.1[[l]],
#         res.2[[l]]
#       )
#     }
#   )
#   names(res) <- names(res.1)
#   res
# }

# rbind_results <- function(...) {
#   Reduce(rbind_2results,list(...))
#
# }

## aggregating results from RDS files ----
# aggregateRDS <- function(folder.path,include.dists = FALSE) {
#   list.files(folder.path,full.names = TRUE) |>
#     lapply(readRDS) %>%
#     {if (include.dists) . else lapply(.,"[",j = -c("dists"))} |>
#     rbindlist()
# }

## Extracting results ----
extract_column <- function(dt,colname) {
  dt |>
    pull(colname) |>
    rbindlist()
}

## Plotting ----
### tools ----
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

### graphs ----
#### Distances & divergences ----
plot_distance <- function(data,dist,ylab) {
  ggplot(data = data,aes_string(x = "type",y = dist,colour = "type",fill = "type"))+
    geom_point(aes_string(group = "interaction(i,j,type)"),alpha = 0.2,position = position_jitterdodge())+
    geom_boxplot(aes_string(group = "interaction(i,j,type)"),alpha = 0.4)+
    geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
    geom_vline(xintercept = 1.5,size = 1.2)+
    scale_x_discrete(labels = x.labs)+
    scale_fill_manual(values = x.fill)+
    scale_colour_manual(values = x.fill)+
    xlab("")+ylab(label = ylab)+
    guides(colour = "none",fill = "none")+
    coord_flip()+
    theme_minimal(12)+theme(axis.text.y = element_text(face = x.labs.face))
}

plot_distance_distrib <- function(data,dist,xlab,x.lims) {
  ggplot(data = data,aes_string(x = dist,y = "type",colour = "type",fill = "type"))+
    geom_density_ridges(aes_string(y = "type",group = "interaction(i,j,type,group.number)"),
                        fill = NA,scale = .8)+
    geom_density_ridges(aes_string(y = "type",group = "interaction(i,j,type,group.number)"),
                        colour = NA,alpha = 0.005,scale = .8)+
    geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
    geom_hline(yintercept = 1.9,size = 1.2)+
    scale_y_discrete(labels = x.labs)+
    scale_x_continuous(limits = x.lims)+
    scale_fill_manual(values = x.fill)+
    scale_colour_manual(values = paste0(x.fill,"25"))+
    ylab("")+xlab(label = xlab)+
    guides(colour = "none",fill = "none")+
    theme_minimal(12)+theme(axis.text.y = element_text(face = x.labs.face))
}


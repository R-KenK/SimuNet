library(data.table)
library(igraph)
library(ggplot2)
library(dplyr)
library(pbapply)
library(ggridges)
devtools::load_all(".")
arrow::set_cpu_count(1)

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
  M[lower.tri(M,diag = TRUE)] <- 0L
  M
}

is_sameClique <- function(group) {
  expand.grid(i = 1:nrow(group),j = 1:nrow(group)) |>
    subset(i < j) |>
    {\(x) {x$same.clique <- group$clique[x$i] == group$clique[x$j];x}}() |> {
      \(x) {
        M <- matrix(NA,nrow(group),nrow(group),
                    dimnames = list(as.character(1:nrow(group)),as.character(1:nrow(group)))
        )
        M[cbind(x$i,x$j)] <- x$same
        M
      }
    }()
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
  M <-
    {within + out} |>
    draw_matbinom(samp.eff = samp.eff)
  mode(M) <- "integer"
  M
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
    param.list[,by = .(n,netgen_name,n.rep,group.number),
               `:=`(netgen_args = list(
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
    # generate one other group per group.number & group.rep
    {
      \(.) .[,by = .(n,netgen_name,n.rep,group.number),
             `:=`(
               netgen_other = list(
                 group = data.table(
                   name = as.character(1:n),
                   value = runif(n),
                   rank = sample(1:n,n),
                   gregariousness = runif(n),
                   clique = sample(LETTERS[1:(rpois(1,1) + 1)],n,replace = TRUE)
                 )
               )
             )
      ]
    }() |>
    {
      \(.) .[,by = .(n,netgen_name,samp.eff,n.rep,group.number,group.rep),
             `:=`(
               fixed.rand.prob = list(runif(n = (n^2 - n) / 2))
             )
      ]
    }() |>
    {
      \(.) .[
        ,
        `:=`(netgen_args = mapply(list,samp.eff = samp.eff,group = netgen_args,
                                  fixed.rand.prob = fixed.rand.prob,SIMPLIFY = FALSE),
             netgen_other = mapply(list,samp.eff = samp.eff,group = netgen_other,SIMPLIFY = FALSE)
        )
      ]
    }() |>
    # aggregates into a list the arguments to be passed to the netgen_fun
    {\(.) .[,fixed.rand.prob := NULL]}()

  param.list
}

## Network data extraction: WIP ----
compute.EV <- function(M) {
  M |>
    igraph::graph.adjacency(mode = "upper",weighted = TRUE)  %>%
    {igraph::eigen_centrality(.)$vector}
}

## Misc. ----
pastetmp <- function(path) {
  if (substr(path,nchar(path),nchar(path)) == "/")
    paste0(substr(path,1,nchar(path) - 1),"tmp")
  else
    paste0(path,"tmp")
}

set.seed.alpha <- function(x) {
  # from https://stackoverflow.com/a/10913336/11939996
  hexval <- paste0("0x",digest::digest(x,"crc32"))
  intval <- type.convert(hexval,as.is = TRUE) %% .Machine$integer.max
  intval
}

## Simulating weighted adjacency matrices ----
run_simulation_single <- function(r,
                                  edgeDT.path = ".WIP/simulation.data/edgeDT/") {
  param.list.row  <- param.list[r,]
  n               <- param.list$n[r];
  samp.eff        <- param.list$samp.eff[r];
  n.rep           <- param.list$n.rep[r];
  names           <- param.list$netgen_args[[r]]$group$name;
  fixed.rand.prob <- param.list$netgen_args[[r]]$fixed.rand.prob;
  f               <- param.list$netgen_fun[[r]];
  args            <- param.list$netgen_args[[r]];
  other           <- param.list$netgen_other[[r]];
  seed            <- param.list$seed[[r]];

  set.seed(seed)

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
      SimuNet::simunet(samp.effort = samp.eff,mode = "upper",n.scans = samp.eff,
                       alpha.prior = 1,beta.prior = 1) |>
      SimuNet::sum_scans()
  )
  param.list.row[
    ,netgen_output := list(
      {
        netgen_output <-
          data.table::data.table(
            type = factor(
              c(
                "real",
                "real.bis",
                "SimuNet",
                "other",
                "ER",
                "fixed.rand",
                "total.rand"
              ),
              levels =
                c(
                  "real",
                  "real.bis",
                  "SimuNet",
                  "other",
                  "ER",
                  "fixed.rand",
                  "total.rand"
                )
            ),
            dist = list(
              real,
              real.bis,
              SimuNet,
              other,
              ER,
              fixed.rand,
              total.rand
            )
          )
        netgen_output[,dist := lapply(dist,adjacencies_to_dt)] |>
          tidyfast::dt_unnest(dist)
      }
    )
  ] |>
    tidyfast::dt_unnest(netgen_output) |>
    arrow::write_dataset(path = edgeDT.path.tmp,
                         partitioning = c("netgen_name","n","samp.eff","group.number","group.rep")
    )
  message("Simulations done!")
}

run_simulations <- function(param.list,n.cores = 7,
                            edgeDT.path = ".WIP/simulation.data/edgeDT/",
                            edgeDT.path.tmp = pastetmp(edgeDT.path),
                            delete.tmp = TRUE,
                            partitioning.vec = c("netgen_name","n","samp.eff")) {
  message("Generating replicable RNG seeds...")
  param.list[
    ,
    seed :=
      paste(netgen_name,n,samp.eff,group.number,group.rep,sep = ".") |>
      sapply(set.seed.alpha)
  ]
  message("Creating parallel workers...")
  cl <- parallel::makeCluster(n.cores)
  on.exit({parallel::stopCluster(cl);rm(cl);gc()})
  parallel::clusterExport(
    cl,
    list(
      "param.list",
      "run_simulation_single",
      "adjacencies_to_dt",
      "edgeDT.path.tmp",
      "is_sameClique","draw_matbinom","greg_product",
      "ER_generate_asso",
      "fixedRandom_generate_asso",
      "totalRandom_generate_asso"
    ),
    envir = environment()
  )
  message("Running the simulations...")
  pbapply::pblapply(
    X = 1:nrow(param.list),
    FUN = run_simulation_single,
    edgeDT.path = edgeDT.path.tmp,
    cl = cl
  )
  message("Flushing parallel workers...")
  parallel::stopCluster(cl);rm(cl);gc();on.exit()
  query_edgeDT(edgeDT.path.tmp) |>
    arrow::write_dataset(path = edgeDT.path,
                         partitioning = partitioning.vec)
  if (delete.tmp)
    unlink(edgeDistanceDT.path.tmp,recursive = TRUE)

  message("Simulations done!")
}

## Calculating edge weight distribution distances ----
### formatting edge weight data to measure distances ----
measure_distances_single <-
  function(r,
           edgeDT.path = ".WIP/simulation.data/edgeDT/",
           edgeDistanceDT.path = ".WIP/simulation.data/edgeDistanceDT/") {
    .netgen_name  <- dist.param$netgen_name[r]
    .n            <- dist.param$n[r]
    .samp.eff     <- dist.param$samp.eff[r]
    .group.number <- dist.param$group.number[r]
    .group.rep    <- dist.param$group.rep[[r]]
    dt <-
      arrow::open_dataset(".WIP/simulation.data/edgeDT/") |>
      dplyr::select(-n.rep,-rep) |>
      dplyr::filter(netgen_name %in% .netgen_name &
                      n %in% .n &
                      samp.eff %in% .samp.eff &
                      group.number %in% .group.number &
                      group.rep %in% .group.rep
                      ) |>
      dplyr::collect() |>
      tidyfast::dt_nest(netgen_name,n,samp.eff,group.number,group.rep,type,i,j,.key = "weight") |>
      {
        \(dt) dt[
          ,possible := lapply(samp.eff,\(s) data.table::data.table(possible = 0:s))
        ][
          ,weight.mat :=
            mapply(\(w,p) factor(w$weight,levels = p$possible),
                   w = weight,p = possible,SIMPLIFY = FALSE) |>
            lapply(table) |>
            {\(.) mapply(\(p,w) cbind(p,weight.prob = c(w / sum(w))),
                         p = possible,w = .,SIMPLIFY = FALSE)}()
        ][
          ,weight.mean := lapply(weight,"[[","weight") |> sapply(mean)
        ][
          ,`:=`(weight =  NULL,possible = NULL)
        ][]
      }() |>
      tidyfast::dt_unnest(weight.mat)

    dt[type == "real"][,type := NULL][] |>
      data.table::setnames(
        old = c("weight.mean","weight.prob"),
        new = c("real.weight.mean","real.weight.prob")) |>
      data.table::merge.data.table(
        x = dt,
        by = c("netgen_name","n","samp.eff","group.number","group.rep","i","j","possible")) |>
      subset(type != "real") |>
      {\(.) {
        .[order(netgen_name,n,samp.eff,group.number,group.rep,type,i,j,possible)
        ][
          ,mean.diff := abs(weight.mean - real.weight.mean)
        ][
          ,c("KS.stat","KS.p") :=
            suppressWarnings(
              ks.test(weight.prob,real.weight.prob,exact = FALSE)[c("statistic","p.value")]
            )
          ,by = .(netgen_name,n,samp.eff,group.number,group.rep,type,i,j)
        ][
          ,.(
            KL    = suppressMessages(philentropy::KL(rbind(weight.prob,real.weight.prob))),
            JS    = suppressMessages(philentropy::JSD(rbind(weight.prob,real.weight.prob)))#,
            # EMD.e = emdist::emd(cbind(possible,weight.prob),
            #                     cbind(possible,real.weight.prob),dist = "euclidean"),
            # EMD.m = emdist::emd(cbind(possible,weight.prob),
            #                     cbind(possible,real.weight.prob),dist = "manhattan")
          )
          ,by = .(netgen_name,n,samp.eff,group.number,group.rep,type,i,j,weight.mean,KS.stat,KS.p)
        ][
          ,reference := "real"
        ]
      }}() |>
      data.table::setcolorder(c("netgen_name","n","samp.eff",
                                "group.number","group.rep",
                                "reference","type","i","j")) |>
      arrow::write_dataset(path = edgeDistanceDT.path,
                           partitioning = c("netgen_name","n","samp.eff","group.number","group.rep"))
    rm(dt);gc()
    message("distances measured!")
  }

measure_distances <- function(param.list,
                                  edgeDT.path = ".WIP/simulation.data/edgeDT/",
                                  edgeDistanceDT.path = ".WIP/simulation.data/edgeDistanceDT/",
                                  edgeDistanceDT.path.tmp = pastetmp(edgeDistanceDT.path),
                                  delete.tmp = TRUE,
                                  partitioning.vec = c("netgen_name","n","samp.eff"),
                                  n.each,n.chunks = 3L,n.cores = 7L) {
  n.each <- max(param.list$group.rep)
  dist.param <-
    param.list |>
    dplyr::select(n,netgen_name,samp.eff,group.number) |>
    dplyr::group_by(n,netgen_name,samp.eff,group.number) |>
    {\(.) suppressMessages(dplyr::summarise(.))}() |>
    data.table::setDT() |>
    {\(.) lapply(parallel::splitIndices(n.each,n.chunks),\(v) {.$group.rep <- list(v);.})}() |>
    data.table::rbindlist() |>
    {\(.) .[order(netgen_name,n,samp.eff,group.number)]}()

  message("Creating parallel workers...")
  cl <- parallel::makeCluster(n.cores)
  on.exit({parallel::stopCluster(cl);rm(cl);gc()})
  parallel::clusterExport(cl,list("dist.param","measure_distances_single"),envir = environment())

  message("Measuring distances...")
  pbapply::pblapply(
    1:nrow(dist.param),
    measure_distances_single,
    edgeDT.path = edgeDT.path,
    edgeDistanceDT.path = edgeDistanceDT.path.tmp,
    cl = cl
  )
  message("Flushing parallel workers...")
  parallel::stopCluster(cl);rm(cl);gc();on.exit()

  message("Aggregating individual .parquet files...")
  query_edgeDistanceDT(edgeDistanceDT.path.tmp) |>
    arrow::write_dataset(path = edgeDistanceDT.path,
                         partitioning = partitioning.vec)

  # param.list.aggregated <-
  #   param.list |>
  #   dplyr::select(-group.number,-group.rep) |>
  #   dplyr::group_by(n,netgen_name,samp.eff) |>
  #   {\(.) suppressMessages(dplyr::summarise(.))}() |>
  #   data.table::setDT()
  # aggregate_parquets(
  #   param.list.aggregated = param.list.aggregated,
  #   dataset = "edgeDistanceDT",
  #   sources.path = edgeDistanceDT.path.tmp,
  #   new.path = edgeDistanceDT.path,
  #   partitioning.vec = partitioning.vec,
  #   n.cores = n.cores
  # )
  if (delete.tmp)
    unlink(edgeDistanceDT.path.tmp,recursive = TRUE)
  message("Distances measured!")
}

### measure edge weight distributions distances ----
# measure_distances_single <-
#   function(r,
#            edgeDistanceData.path = ".WIP/simulation.data/edgeDistanceData/",
#            edgeDistanceDT.path = ".WIP/simulation.data/edgeDistanceDT/") {
#     .netgen_name  <- dist.param$netgen_name[r]
#     .n            <- dist.param$n[r]
#     .samp.eff     <- dist.param$samp.eff[r]
#     .group.number <- dist.param$group.number[r]
#     .group.rep    <- dist.param$group.rep[[r]]
#     arrow::open_dataset(sources = edgeDistanceData.path,format = "parquet") |>
#       dplyr::filter(netgen_name %in% .netgen_name &
#                       n %in% .n &
#                       samp.eff %in% .samp.eff &
#                       group.number %in% .group.number &
#                       group.rep %in% .group.rep
#                     ) |>
#       dplyr::collect() |>
#       {\(.) {
#         .[order(netgen_name,n,samp.eff,group.number,group.rep,type,i,j,possible)
#         ][
#           ,mean.diff := abs(weight.mean - real.weight.mean)
#         ][
#           ,c("KS.stat","KS.p") :=
#             suppressWarnings(
#               ks.test(weight.prob,real.weight.prob,exact = FALSE)[c("statistic","p.value")]
#             )
#           ,by = .(netgen_name,n,samp.eff,group.number,group.rep,type,i,j)
#         ][
#           ,.(
#             KL    = suppressMessages(philentropy::KL(rbind(weight.prob,real.weight.prob))),
#             JS    = suppressMessages(philentropy::JSD(rbind(weight.prob,real.weight.prob)))#,
#             # EMD.e = emdist::emd(cbind(possible,weight.prob),
#             #                     cbind(possible,real.weight.prob),dist = "euclidean"),
#             # EMD.m = emdist::emd(cbind(possible,weight.prob),
#             #                     cbind(possible,real.weight.prob),dist = "manhattan")
#           )
#           ,by = .(netgen_name,n,samp.eff,group.number,group.rep,type,i,j,weight.mean,KS.stat,KS.p)
#         ][
#           ,reference := "real"
#         ]
#       }}() |>
#       data.table::setcolorder(c("netgen_name","n","samp.eff",
#                     "group.number","group.rep",
#                     "reference","type","i","j")) |>
#       arrow::write_dataset(path = edgeDistanceDT.path,
#                            partitioning = c("netgen_name","n","samp.eff","group.number","group.rep"))
#     rm(dt);gc()
#     NULL
#   }
#
# measure_distances <- function(param.list,
#                               edgeDistanceData.path = ".WIP/simulation.data/edgeDistanceData/",
#                               edgeDistanceDT.path = ".WIP/simulation.data/edgeDistanceDT/",
#                               edgeDistanceDT.path.tmp = pastetmp(edgeDistanceDT.path),
#                               delete.tmp = TRUE,
#                               partitioning.vec = c("netgen_name","n","samp.eff"),
#                               n.each,n.chunks = 3L,n.cores = 7L) {
#   dist.param <-
#     param.list |>
#     dplyr::select(n,netgen_name,samp.eff,group.number) |>
#     dplyr::group_by(n,netgen_name,samp.eff,group.number) |>
#     {\(.) suppressMessages(dplyr::summarise(.))}() |>
#     data.table::setDT() |>
#     {\(.) lapply(parallel::splitIndices(n.each,n.chunks),\(v) {.$group.rep <- list(v);.})}() |>
#     data.table::rbindlist() |>
#     {\(.) .[order(netgen_name,n,samp.eff,group.number)]}()
#
#   message("Creating parallel workers...")
#   cl <- parallel::makeCluster(n.cores)
#   on.exit({parallel::stopCluster(cl);rm(cl);gc()})
#   parallel::clusterExport(cl,list("dist.param","measure_distances_single"),envir = environment())
#
#   message("Measuring distances...")
#   pbapply::pblapply(
#     1:nrow(dist.param),
#     measure_distances_single,
#     edgeDistanceData.path = edgeDistanceData.path,
#     edgeDistanceDT.path = edgeDistanceDT.path.tmp,
#     cl = cl
#   )
#   message("Flushing parallel workers...")
#   parallel::stopCluster(cl);rm(cl);gc();on.exit()
#
#   query_edgeDistanceDT(edgeDistanceDT.path.tmp) |>
#     arrow::write_dataset(path = edgeDistanceDT.path,
#                          partitioning = partitioning.vec)
#   if (delete.tmp)
#     unlink(edgeDistanceDT.path.tmp,recursive = TRUE)
#   message("Distances measured!")
# }

## Aggregating edge weight data ----
adjacencies_to_dt <- function(adj.array) {
  .n     <- dim(adj.array)[1]
  .n.rep <- dim(adj.array)[3]

  expand.grid(i = 1:.n,j = 1:.n,rep = 1:.n.rep) |>
    data.table::setDT() |>
    cbind(weight = c(adj.array)) |>
    subset(i < j)
}

## Aggregate parquet files ----
aggregate_parquets_single <- function(r,
                                      param.list.aggregated = param.list.aggregated,
                                      dataset = c("edgeDT","edgeDistanceDT","edgeDistanceData"),
                                      sources.path,
                                      new.path,
                                      partitioning.vec = c("netgen_name","n","samp.eff")) {
  dataset <- match.arg(dataset)
  query_fun <- switch(dataset,"edgeDT" = query_edgeDT,
                      "edgeDistanceDT" = query_edgeDistanceDT,
                      "edgeDistanceData" = query_edgeDistanceData
  )
  .n             <- param.list.aggregated$n[r]
  .samp.eff      <- param.list.aggregated$samp.eff[r]
  # .group.number  <- param.list.aggregated$group.number[r]

  query_fun(sources.path) |>
    dplyr::filter(n == .n & samp.eff == .samp.eff) |># & group.number == .group.number) |>
    arrow::write_dataset(path = new.path,
                         partitioning = partitioning.vec)
  message("Aggregation done!")
  "Aggregation done!"
}

aggregate_parquets <- function(param.list.aggregated,
                               dataset = c("edgeDT","edgeDistanceDT","edgeDistanceData"),
                               sources.path,
                               new.path,
                               partitioning.vec = c("netgen_name","n","samp.eff"),
                               n.cores = 7L) {
  message("Creating parallel workers...")
  cl <- parallel::makeCluster(n.cores)
  on.exit({parallel::stopCluster(cl);rm(cl);gc()})
  parallel::clusterExport(cl,
                          list("param.list.aggregated",
                               "aggregate_parquets_single",
                               "query_edgeDT",
                               "query_edgeDistanceDT",
                               "query_edgeDistanceData"),
                          envir = environment())
  message("Aggregating parquet files...")
  pbapply::pblapply(
    X = 1:nrow(param.list.aggregated),
    FUN = aggregate_parquets_single,
    param.list.aggregated = param.list.aggregated,
    dataset = dataset,
    sources.path = sources.path,
    new.path = new.path,
    partitioning.vec = partitioning.vec,
    cl = cl
  )
}


## Retrieve from Arrow's datasets ----
complete_edgedt <- function(dt,.group.rep,n,n.rep) {
  expand.grid(group.rep = .group.rep,rep = 1:n.rep,i = 1:n,j = 1:n,weight = 0L) |>
    setDT() |>
    subset(i >= j) |>
    rbind(dt) |>
    arrange(rep,group.rep,j,i)
}

reconstruct_adjacencies <- function(.netgen_name,
                                    .n,
                                    .samp.eff,
                                    .type,
                                    .group.number,
                                    .group.rep,
                                    n.rep = 105L,
                                    edgeDT.path = ".WIP/simulation.data/edgeDT/") {
  query_edgeDT(edgeDT.path = edgeDT.path) |>
    dplyr::filter(
      netgen_name  %in% .netgen_name,
      n            %in% .n,
      samp.eff     %in% .samp.eff,
      type         %in% .type,
      group.number %in% .group.number,
      group.rep    %in% .group.rep) |>
    dplyr::select(group.rep,rep,i,j,weight) |>
    collect() |>
    complete_edgedt(.group.rep = .group.rep,n = .n,n.rep = n.rep * length(.group.rep)) |>
    dplyr::pull(weight) |>
    array(c(.n,.n,n.rep * length(.group.rep)),dimnames = list(as.character(1:.n),as.character(1:.n),NULL))
}

query_edgeDT <- function(edgeDT.path = ".WIP/simulation.data/edgeDT/") {
  arrow::open_dataset(sources = edgeDT.path) |>
    dplyr::relocate(c("netgen_name","n","samp.eff",
                      "group.number","group.rep",
                      "type","n.rep","rep","i","j","weight")) |>
    dplyr::arrange(n,samp.eff,group.number,group.rep,n.rep,rep,i,j)
}

query_edgeDistanceDT <- function(edgeDistanceDT.path = ".WIP/simulation.data/edgeDistanceDT/") {
  arrow::open_dataset(sources = edgeDistanceDT.path) |>
    dplyr::relocate(c("netgen_name","n","samp.eff",
                      "group.number","group.rep",
                      "reference","type","i","j",
                      "weight.mean","KS.stat","KS.p",
                      "KL","JS")) |>
    dplyr::arrange(n,samp.eff,group.number,group.rep,i,j)
}

query_edgeDistanceData <- function(edgeDistanceData.path = ".WIP/simulation.data/edgeDistanceData/") {
  arrow::open_dataset(sources = edgeDistanceData.path) |>
    dplyr::relocate(c("netgen_name","n","samp.eff",
                      "group.number","group.rep",
                      "type","i","j",
                      "weight.mean","weight.prob",
                      "real.weight.mean","weight.prob")) |>
    dplyr::arrange(n,samp.eff,group.number,group.rep,i,j)
}

## Plotting ----
### tools ----
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

### graphs ----
#### Distances & divergences ----
x.labs <- c(
  "...the same network (expected null distance)",
  paste0("... ",c("the SimuNet's approximation",
                  "another network\n(similar generation, different group)",
                  "an Erdős–Rényi graph (p = 0.5)",
                  "a Random network (fixed prob.)",
                  "a Random network (variable prob.)"))
)

x.labs.face <- c(
  "bold",
  rep("plain",length(x.labs) - 1)
)

x.fill <- c("#FF6347","#4169E1",gg_color_hue(5)[2:5]) # "tomato" and "royalblue"

plot_distance <- function(data,dist,ylab,geom = c("boxplot","density"),
                          x.lims = c(NA,NA),.group = "interaction(i,j,type,group.number)",.alpha = 0.005) {
  geom <- match.arg(geom)
  switch(
    geom,
    "boxplot" =
      ggplot(data = data,aes_string(x = "type",y = dist,colour = "type",fill = "type"))+
      geom_point(aes_string(group = .group),alpha = 0.2,position = position_jitterdodge())+
      geom_boxplot(aes_string(group = .group),alpha = 0.4)+
      geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
      geom_vline(xintercept = 1.5,size = 1.2)+
      scale_x_discrete(labels = x.labs)+
      scale_fill_manual(values = x.fill)+
      scale_colour_manual(values = x.fill)+
      xlab("")+ylab(label = ylab)+
      guides(colour = "none",fill = "none")+
      coord_flip()+
      theme_minimal(12)+theme(axis.text.y = element_text(face = x.labs.face)),
    "density" =
      ggplot(data = data,aes_string(x = dist,y = "type",colour = "type",fill = "type"))+
      geom_density_ridges(aes_string(y = "type",group = .group),
                          fill = NA,scale = .8)+
      geom_density_ridges(aes_string(y = "type",group = .group),
                          colour = NA,alpha = .alpha,scale = .8)+
      geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
      geom_hline(yintercept = 1.9,size = 1.2)+
      scale_y_discrete(labels = x.labs)+
      scale_x_continuous(limits = x.lims)+
      scale_fill_manual(values = x.fill)+
      scale_colour_manual(values = paste0(x.fill,"25"))+
      ylab("")+xlab(label = ylab)+
      guides(colour = "none",fill = "none")+
      theme_minimal(12)+theme(axis.text.y = element_text(face = x.labs.face))
  )
}

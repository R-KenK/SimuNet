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
run_simulation_single <- function(r) {
  n               <- param.list$n[r];
  samp.eff        <- param.list$samp.eff[r];
  n.rep           <- param.list$n.rep[r];
  names           <- param.list$netgen_args[[r]]$group$name;
  fixed.rand.prob <- param.list$netgen_args[[r]]$fixed.rand.prob;
  f               <- param.list$netgen_fun[[r]];
  args            <- param.list$netgen_args[[r]];
  other           <- param.list$netgen_other[[r]];

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
  data.table::data.table(
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

run_simulations <- function(param.list,n.cores = 7) {
  message("Creating parallel workers...")
  cl <- parallel::makeCluster(n.cores)
  on.exit({parallel::stopCluster(cl);rm(cl);gc()})
  parallel::clusterExport(
    cl,
    list(
      "param.list",
      "run_simulation_single",
      "is_sameClique","draw_matbinom","greg_product",
      "ER_generate_asso",
      "fixedRandom_generate_asso",
      "totalRandom_generate_asso"
    )
  )
  message("Running the simulations...")
  param.list[
    ,
    netgen_output :=
      pbapply::pblapply(
        X = 1:nrow(param.list),
        FUN = run_simulation_single,
        cl = cl
      )
  ]
  param.list[]
}

agregate_edgeDT_single <- function(r,dataset.path = ".WIP/simulation.data/edgeDT/") {
  n            <- param.list$n[r]
  samp.eff     <- param.list$samp.eff[r]
  netgen_name  <- param.list$netgen_name [r]
  n.rep        <- param.list$n.rep[r]
  group.number <- param.list$group.number[r]
  group.rep    <- param.list$group.rep[r]
  real         <- param.list$netgen_output[[r]][type == "real"]$dist[[1]][[1]]
  real.bis     <- param.list$netgen_output[[r]][type == "real.bis"]$dist[[1]][[1]]
  other        <- param.list$netgen_output[[r]][type == "other"]$dist[[1]][[1]]
  ER           <- param.list$netgen_output[[r]][type == "ER"]$dist[[1]][[1]]
  fixed.rand   <- param.list$netgen_output[[r]][type == "fixed.rand"]$dist[[1]][[1]]
  total.rand   <- param.list$netgen_output[[r]][type == "total.rand"]$dist[[1]][[1]]
  SimuNet      <- param.list$netgen_output[[r]][type == "SimuNet"]$dist[[1]][[1]]
  edge.dt <-
    expand.grid(i = 1:n,j = 1:n,rep = 1:n.rep) |>
    subset(i < j) |>
    {
      \(x) rbind(
        cbind(x,n = n,samp.eff = samp.eff,n.rep = n.rep,type = "real",
              weight = real[as.matrix(x)]),
        cbind(x,n = n,samp.eff = samp.eff,n.rep = n.rep,type = "real.bis",
              weight = real.bis[as.matrix(x)]),
        cbind(x,n = n,samp.eff = samp.eff,n.rep = n.rep,type = "other",
              weight = other[as.matrix(x)]),
        cbind(x,n = n,samp.eff = samp.eff,n.rep = n.rep,type = "ER",
              weight = ER[as.matrix(x)]),
        cbind(x,n = n,samp.eff = samp.eff,n.rep = n.rep,type = "fixed.rand",
              weight = fixed.rand[as.matrix(x)]),
        cbind(x,n = n,samp.eff = samp.eff,n.rep = n.rep,type = "total.rand",
              weight = total.rand[as.matrix(x)]),
        cbind(x,n = n,samp.eff = samp.eff,n.rep = n.rep,type = "SimuNet",
              weight = SimuNet[as.matrix(x)])
      )
    }() |>
    data.table()
  edge.dt$type <-
    edge.dt$type |>
    factor(levels = c("real","real.bis","SimuNet","other","ER","fixed.rand","total.rand"))
  cbind(edge.dt,netgen_name = netgen_name,group.number = group.number,group.rep = group.rep) |>
    arrow::write_dataset(path = dataset.path,
                         partitioning = c("netgen_name","n","samp.eff","group.number","group.rep")
    )
  NULL
}

aggregate_edgeDT <- function(param.list,n.cores = 7) {
  message("Creating parallel workers...")
  cl <- parallel::makeCluster(n.cores)
  on.exit({parallel::stopCluster(cl);rm(cl);gc()})
  parallel::clusterExport(
    cl,
    list(
      "param.list",
      "agregate_edgeDT_single",
      "data.table"
    )
  )
  message("Aggregating edge weights into Arrow dataset...")
  pbapply::pblapply(
    X = 1:nrow(param.list),
    FUN = agregate_edgeDT_single,
    cl = cl
  )
}

calculate_edgeDistanceDT <- function(param.list,cl) {
  param.list <- prepare_for_distances(param.list)
  measure_all_distances(param.list,cl = cl)
}

## Aggregating experiment results ----

# ## Extracting results ----
# extract_column <- function(dt,colname) {
#   dt |>
#     pull(colname) |>
#     rbindlist()
# }
#
# get_from_RDS <- function(path,dt.name) {
#   dt <- readRDS(path)
#   dt[[dt.name]] |>
#     data.table::rbindlist()
# }


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


# Load packages and custom functions ----
source(".WIP/experiments_on_modelling_networks_generated_by_complex_processes_tools.R")

# Case study ----
## Generating main variables ----
set.seed(42)

## Simulations ----
n.rep <- 105L
n.group <- 3L
n.each <- 49L

### Running the simulations ----
param.n <-
  #seq(6,30,by = 3)
  c(5,8,10,15,20)
param.samp.eff <-
  # seq(20,500,by = 20)
  c(10,25,50,75,100,150,250,500)

param.netgen <-
  list(
    # data.table(netgen_name = "HRCMV",netgen_fun = list(highRankChooseMostValuable_generate_asso)),
    data.table(netgen_name = "GWC",netgen_fun = list(gregWithinClique_generate_asso))
  )

param.list <-
  generate_paramList(param.n = param.n,
                     param.samp.eff = param.samp.eff,
                     param.netgen = param.netgen,
                     n.rep = n.rep,
                     n.group = n.group,
                     n.each = n.each
  )[group.number <= 2 & group.rep <= 28]
param.list[]

### Running simulations ----
start.time <- Sys.time()
start.time
param.list <- run_simulations(param.list,n.cores = 7)
end.time <- Sys.time()
end.time
end.time - start.time
param.list

aggregate_edgeDT(param.list,n.cores = 7)

arrow::open_dataset(".WIP/simulation.data/edgeDTtest/") |>
  pull(group.rep) |> unique() |> sort()

arrow::open_dataset(".WIP/simulation.data/edgeDT/") |>
  group_by(n,samp.eff,group.number,group.rep) |>
  collect() |>
  count(name = "nrow") |>
  group_by(n,nrow) |>
  count(name = "m")

arrow::open_dataset(".WIP/simulation.data/edgeDTtest/") |>
  group_by(n,samp.eff,group.number,group.rep) |>
  collect() |>
  count(name = "nrow") |>
  group_by(n,nrow) |>
  count(name = "m")

results <- prepare_for_distances(results)

snow::clusterCall(cl,\() {rm(results, pos = globalenv()); gc()})
snow::clusterExport(cl,list("results"),envir = environment())

results <- measure_all_distances(results,cl = cl)

end.time <- Sys.time()
end.time
end.time - start.time
snow::stopCluster(cl);rm(cl)
results

# saveRDS(results,".WIP/simulation.data/results.GWC.n_8.10.15.seff_75.100.150.rds")

# Results wrangling ----
## Aggregate RDS files ----
filepaths <- list.files(".WIP/simulation.data/",pattern = "GWC",full.names = TRUE)

## Aggregating edge.distances.dt ----
cl <- snow::makeCluster(7)
snow::clusterExport(cl,list("filepaths","get_from_RDS"))
edge.distances.dt <- pbapply::pblapply(filepaths,get_from_RDS,dt.name = "edge.distances.dt",cl = cl) |>
  rbindlist(use.names = TRUE)
snow::stopCluster(cl);rm(cl)

## Saving in optimal file type ----
edge.distances.dt |>
  arrow::write_dataset(
    path = ".WIP/simulation.data/edgeDistanceDT/",
    format = "parquet",
    partitioning = c("n")
  )

## Retrieving from optimal file type ----
col.labs <- c(
  "Itself (expected null distance)",
  "SimuNet",
  "another network (similar generation)",
  "an Erdős–Rényi graph (p = 0.5)",
  "a Random network (fixed prob.)",
  "a Random network (variable prob.)"
)

arrow::open_dataset(sources = ".WIP/simulation.data/edgeDistanceDT/") |>
  select(n, samp.eff, group.number, reference, type, i, j, KL) |>
  group_by(n, samp.eff, group.number, reference, type) |>
  summarise(
    KL = mean(KL),
    inf   = quantile(KL, probs = 0.25),
    sup   = quantile(KL, probs = 0.75)
  ) |>
  collect() |>
  subset(reference == "real") |>
  mutate(linet = ifelse(type %in% c("real.bis","SimuNet"),"dashed",NA)) |>
  mutate(alph = ifelse(type %in% c("real.bis","SimuNet"),"1","2")) |>
  ggplot(aes(samp.eff,KL,colour = type,fill = type))+
  facet_grid(paste0("repeatition n°",group.number) ~ paste0("n = ",n))+
  geom_ribbon(aes(ymin = inf,ymax = sup,lty = linet,alpha = alph))+
  geom_line()+
  geom_point()+
  scale_fill_manual(labels = col.labs,values = c("tomato","royalblue",gg_color_hue(7)[2:5]))+
  scale_colour_manual(labels = col.labs,values = c("tomato","royalblue",gg_color_hue(7)[2:5]))+
  scale_linetype_manual(values = c("dotted"))+
  scale_alpha_manual(values = c(.3,.1))+
  geom_vline(xintercept = 0,colour = "black")+
  geom_hline(yintercept = 0,colour = "black")+
  scale_x_continuous(expand = c(0,NA))+
  scale_y_continuous(expand = c(0,NA))+
  guides(linetype = "none",alpha = "none")+
  xlab("Sampling effort")+ylab("Kullback-Leibler divergence")+
  labs(colour = "Real network compared to...",fill = "Real network compared to...")+
  theme_minimal(14)+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

## Aggregating edge.dt: too big in a single data.table... ----
source(".WIP/experiments_on_modelling_networks_generated_by_complex_processes_tools.R")
filepaths <-
  list.files(".WIP/simulation.data/",pattern = "GWC",full.names = TRUE) |>
  grep(pattern = "20_",value = TRUE)

cl <- snow::makeCluster(min(7,length(filepaths)))
snow::clusterExport(cl,list("filepaths","get_from_RDS"))
edge.dt <-
  pbapply::pblapply(filepaths,get_from_RDS,dt.name = "edge.dt",cl = cl) |>
  rbindlist(use.names = TRUE) |>
  subset(n == 20)
snow::stopCluster(cl);rm(cl);gc()

## Saving in optimal file type ----
arrow::write_parquet(edge.dt,sink = ".WIP/simulation.data/parquet/edge.dt.n_20.parquet")
rm(edge.dt);gc()

## Retrieving from optimal file type ----
edge.dt <- arrow::read_parquet(file = ".WIP/simulation.data/test/edge.dt.parquet")

edge.dt <-
  list.files(".WIP/simulation.data/parquet/",pattern = "parquet",full.names = TRUE) |>
  lapply(arrow::read_parquet) |>
  rbindlist()

arrow::write_dataset(edge.dt,
                     path = ".WIP/simulation.data/edgeDT/",format = "parquet",
                     partitioning = c("n","group.number")
)

weidgt.dt <-
  arrow::open_dataset(".WIP/simulation.data/edgeDT") |>
  filter(n %in% c(5,8,10)) |>
  select(i,j,n,samp.eff,type,group.number,group.rep,weight) |>
  group_by(i,j,n,samp.eff,type,group.number,group.rep) |>
  summarise(
    weight = median(weight),
    inf = quantile(weight,probs = 0.025),
    sup = quantile(weight,probs = 0.975)
  ) |>
  dplyr::collect()

setDT(weidgt.dt)

weidgt.dt


n.val <-
  arrow::open_dataset(".WIP/simulation.data/edgeDT") |>
  select(n) |>
  distinct() |>
  collect() |>
  pull()

w.dt <-
  lapply(
    c(5,8,10),
    function(val) {
      arrow::open_dataset(".WIP/simulation.data/edgeDT") |>
        filter(n == val) |>
        select(i,j,n,samp.eff,type,group.number,weight) |>
        mutate(weight = weight / samp.eff) |>
        group_by(i,j,n,samp.eff,type,group.number) |>
        summarise(
          weight = median(weight),
          inf = quantile(weight,probs = 0.025),
          sup = quantile(weight,probs = 0.975)
        ) |>
        dplyr::collect()
    }
  ) |>
  bind_rows()

setDT(w.dt)

type.str <- c("the real network",
              "the real network bis",
              "SimuNet network",
              "another network (similar generation)",
              "an Erdős–Rényi graph (p = 0.5)",
              "a Random network (fixed prob.)",
              "a Random network (variable prob.)")

w.dt |>
  # subset(i <= 5 & j <= 5) |>
  subset(group.number == 1) |>
  subset(n == 10) %>% {
    ggplot(data = .,aes(type,weight, colour = samp.eff, fill = type))+
      geom_line(aes(group = interaction(group.number,samp.eff)),lty = "dashed",alpha = 0.6,position = position_dodge(0.5))+
      geom_errorbar(aes(group = interaction(group.number,samp.eff),ymin = inf,ymax = sup),position = position_dodge(0.5))+
      geom_point(aes(group = interaction(group.number,samp.eff)),position = position_dodge(0.5))+
      facet_grid(i ~ j,scales = "free_y")+
      # geom_density(bw = 0.05,alpha = 0.25,position = "identity",colour = NA)+
      # geom_density(data = .[type %in% c("SimuNet","real")],
      #              bw = 0.05,alpha = 0.15,position = "identity")+
      # geom_vline(data = SimuNet.dt,aes(xintercept = weight),
      #            colour = "tomato",lty = "dashed",size = 1.1)+
      # geom_vline(data = real.dt,aes(xintercept = weight),
      #            colour = "royalblue",lty = "dashed",size = 1.1)+
      # geom_point(data = Adj.obs.dt,aes(x = weight,y = 0),fill = "white",
      #            colour = "tomato",shape = 21,stroke = 1.5,size = 2)+
      # scale_colour_manual(values = c("tomato","royalblue",gg_color_hue(6)[2:6]),labels = type.str)+
      scale_colour_gradient(low = "royalblue",high = "orange")+
      scale_fill_manual(values = c("tomato","royalblue",gg_color_hue(6)[2:6]),labels = type.str)+
      scale_x_discrete(labels = type.str)+
      scale_y_continuous(breaks = c(0,.5,1),labels = c('0','0.5','1'),expand = c(0,NA))+
      labs(colour = "Sampling effort")+
      guides(fill = "none")+
      xlab("Network type")+ylab("Edge weight")+
      geom_hline(yintercept = 0)+
      # geom_vline(xintercept = 0)+
      coord_flip()+
      cowplot::theme_minimal_hgrid(12)+theme()
  }


  arrow::open_dataset(".WIP/simulation.data/edgeDT") |>
  filter(n == 5) |>
  select(i,j,n,samp.eff,type,group.number,group.rep,weight) |>
  mutate(weight = weight / samp.eff) |>
  group_by(i,j,n,samp.eff,type,group.number,group.rep) |>
  summarise(
    weight = median(weight),
    inf = quantile(weight,probs = 0.025),
    sup = quantile(weight,probs = 0.975)
  )

edge.distances.dt <-
  edge.distances.dt[order(n,samp.eff,group.number,group.rep,reference,type,i,j)] %>%
  .[,by = .(n,samp.eff,group.number,reference,type,i,j),
    .(
      mean.diff= median(mean.diff),
      KS.stat  = median(KS.stat),
      KS.p     = median(KS.p),
      KL       = median(KL),
      JS       = median(JS),
      EMD.e    = median(EMD.e),
      EMD.m    = median(EMD.m),
      mean.diff.low= quantile(probs = 0.25,mean.diff),
      KS.stat.low  = quantile(probs = 0.25,KS.stat),
      KS.p.low     = quantile(probs = 0.25,KS.p),
      KL.low       = quantile(probs = 0.25,KL),
      JS.low       = quantile(probs = 0.25,JS),
      EMD.e.low    = quantile(probs = 0.25,EMD.e),
      EMD.m.low    = quantile(probs = 0.25,EMD.m),
      mean.diff.hig= quantile(probs = 0.975,mean.diff),
      KS.stat.hig  = quantile(probs = 0.975,KS.stat),
      KS.p.hig     = quantile(probs = 0.975,KS.p),
      KL.hig       = quantile(probs = 0.975,KL),
      JS.hig       = quantile(probs = 0.975,JS),
      EMD.e.hig    = quantile(probs = 0.975,EMD.e),
      EMD.m.hig    = quantile(probs = 0.975,EMD.m)
    )
  ]

get_from_RDS(".WIP/simulation.data/results.GWC.n_5.8_seff_10.25.50.rds","edge.distance.dt")

test <- readRDS(".WIP/simulation.data/results.GWC.n_5.8_seff_10.25.50.rds")$edge.distances.dt
tast <- test$edge.distance.dt |> rbindlist()
test.2 <- test[,netgen_output := NULL]

results |> object.size() %>% {. / 1024^2}
test$netgen_fun |> as.character() |> as.list() |> as.function()


# Importing results ----
list.files(".WIP/simulation.data/",full.names = TRUE)
data.table::fread(".WIP/simulation.data/results.GWC.n_8.10.15.seff_75.100.150.rds")
results <- readRDS(".WIP/simulation.data/results.GWC.n_8.10.15.seff_75.100.150.rds")
results

# test bench measure distances ----


test.results <- copy(results)
results <- prepare_for_distances(results,"real","SimuNet")
snow::clusterExport(cl,list("results"))
results <- measure_distances_vec(results,cl = cl)
results$distances[[220]]
testSN <- measure_distances2(results,"real","SimuNet",cl = NULL)
testrealB <- measure_distances2(results,"real","real.bis",cl = NULL)
testER <- measure_distances2(results,"real","ER",cl = NULL)

list(testSN,testrealB,testER) |>
  split_seq() |>
  lapply(do.call,what = rbind)

{
  do.call(mapply,list(FUN = rbind,... = .,SIMPLIFY = FALSE))
}

  mapply(FUN = rbind,SIMPLIFY = FALSE)

test.results[,edge.weight := mapply(cbind,
                                  n = n,samp.eff = samp.eff,n.rep = n.rep,group.number = group.number,group.rep = group.rep,edge.weight,SIMPLIFY = FALSE)]
test.results[]
test.results$edge.weight[[220]]

## Plotting distances ----
x.labs <- c(
  "Expected null distance",
  paste0("... ",c("the real network",
                  "another network (similar generation)",
                  "an Erdős–Rényi graph (p = 0.5)",
                  "a Random network (fixed prob.)",
                  "a Random network (variable prob.)"))
)

x.labs.face <- c(
  "bold",
  rep("plain",length(x.labs) - 1)
)

x.fill <- c("#FF6347","#4169E1",gg_color_hue(5)[2:5]) # "tomato" and "royalblue"
results[,edge.distance.dt :=lapply(edge.distance.dt,\(ddt) ddt$type <- ddt$type |>
                                     factor(levels = c("real","real.bis","SimuNet","other","ER","fixed.rand","total.rand"))
)]
results[,type := type |>
          factor(levels = c("real","real.bis","SimuNet","other","ER","fixed.rand","total.rand"))]
### Single case as boxplots ----
KS.stat <-
  results |>
  subset(samp.eff == 75) |>
  subset(group.number == 1) |>
  extract_column("edge.distance.dt") |>
  ggplot(aes_string(x = "type",y = "KS.stat",colour = "type",fill = "type"))+
  geom_point(aes_string(group = "interaction(i,j,type)"),alpha = 0.2,position = position_jitterdodge())+
  geom_boxplot(aes_string(group = "interaction(i,j,type)"),alpha = 0.4)+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_vline(xintercept = 1.5,size = 1.2)+
  # scale_x_discrete(labels = x.labs)+
  scale_fill_manual(values = x.fill)+
  scale_colour_manual(values = x.fill)+
  # xlab("")+ylab(label = ylab)+
  guides(colour = "none",fill = "none")+
  coord_flip()+
  theme_minimal(12)+theme(axis.text.y = element_text(face = x.labs.face))
  plot_distance(dist = "KS.stat",ylab = "Kolmogorov-Smirnoff statistic")

KS.p <-
  results |>
  subset(samp.eff == 75) |>
  subset(group.number == 1) |>
  extract_column("edge.distance.dt") |>
  plot_distance(dist = "KS.p",ylab = "Kolmogorov-Smirnoff p-value")

KL <-
  results |>
  subset(samp.eff == 75) |>
  subset(group.number == 1) |>
  extract_column("edge.distance.dt") |>
  plot_distance(dist = "KL",ylab = "Kullback-Leibler divergence")

JS <-
  results |>
  subset(samp.eff == 75) |>
  subset(group.number == 1) |>
  extract_column("edge.distance.dt") |>
  plot_distance(dist = "JS",ylab = "Jensen-Shannon distance")

EMD.e <-
  results |>
  subset(samp.eff == 75) |>
  subset(group.number == 1) |>
  extract_column("edge.distance.dt") |>
  plot_distance(dist = "EMD.e",ylab = "Earth Mover Distance (Euclidean dist.)")

EMD.m <-
  results |>
  subset(samp.eff == 75) |>
  subset(group.number == 1) |>
  extract_column("edge.distance.dt") |>
  plot_distance(dist = "EMD.m",ylab = "Earth Mover Distance (Manhattan dist.)")

gridExtra::grid.arrange(KS.stat,KS.p,KL,JS,EMD.e,EMD.m,ncol = 2,
                        top = grid::textGrob("SimuNet compared to...",
                                             gp=grid::gpar(fontsize = 16,fontface = "bold")))


### Multiple cases as densities ----
KS.stat <-
  results |>
  subset(n == 10) |>
  subset(samp.eff == 75) |>
  subset(group.number <= 1) |>
  extract_column("edge.distance.dt") |>
  plot_distance_distrib(dist = "KS.stat",xlab = "Kolmogorov-Smirnoff statistic",x.lims = c(0,1))
  # facet_grid(. ~ samp.eff)

KS.p <-
  results |>
  subset(n == 10) |>
  subset(samp.eff == 75) |>
  subset(group.number <= 1) |>
  extract_column("edge.distance.dt") |>
  plot_distance_distrib(dist = "KS.p",xlab = "Kolmogorov-Smirnoff p-value",x.lims = c(0,1))#+
# facet_grid(. ~ samp.eff)

KL <-
  results |>
  subset(n == 10) |>
  subset(samp.eff == 75) |>
  subset(group.number <= 1) |>
  # extract_column("edge.distance.dt") |>
  plot_distance_distrib(dist = "KL",xlab = "Kullback-Leibler divergence",x.lims = c(0,NA))#+
# facet_grid(. ~ samp.eff)

JS <-
  results |>
  subset(n == 10) |>
  subset(samp.eff == 75) |>
  subset(group.number <= 1) |>
  extract_column("edge.distance.dt") |>
  plot_distance_distrib(dist = "JS",xlab = "Jensen-Shannon distance",x.lims = c(0,NA))#+
  # facet_grid(. ~ samp.eff)

EMD.e <-
  results |>
  subset(n == 10) |>
  subset(samp.eff == 75) |>
  subset(group.number <= 1) |>
  extract_column("edge.distance.dt") |>
  plot_distance_distrib(dist = "EMD.e",xlab = "Earth Mover Distance (Euclidean dist.)",x.lims = c(0,NA))#+
# facet_grid(. ~ samp.eff)

EMD.m <-
  results |>
  subset(n == 10) |>
  subset(samp.eff == 75) |>
  subset(group.number <= 1) |>
  extract_column("edge.distance.dt") |>
  plot_distance_distrib(dist = "EMD.m",xlab = "Earth Mover Distance (Manhattan dist.)",x.lims = c(0,NA))#+
# facet_grid(. ~ samp.eff)

gridExtra::grid.arrange(KS.stat,KS.p,KL,JS,EMD.e,EMD.m,ncol = 2,
                        top = grid::textGrob("SimuNet compared to...",
                                             gp=grid::gpar(fontsize = 16,fontface = "bold")))

# Graphic exploration ---------------------------------------------------------------------------------------
## single sim ----
type.str <- c("the real network",
              "the real network bis",
              "SimuNet network",
              "another network (similar generation)",
              "an Erdős–Rényi graph (p = 0.5)",
              "a Random network (fixed prob.)",
              "a Random network (variable prob.)")
n.single <-
  results |>
  subset(samp.eff == 50) |>
  subset(group.number == 1 & group.rep == 1) |>
  pull(n)

Adj.obs.dt <-
  expand.grid(i = 1:n.single,j = 1:n.single) |>
  subset(i < j) |>
  as.matrix() %>%
  {
    cbind(
      .,
      weight = {Adj <- results |>
        subset(samp.eff == 50) |>
        subset(group.number == 1 & group.rep == 1) %>%
        {pull(.,dists)[[1]]} %>%
        {.$real.dist[[1]][,,1]};Adj[.]}
    )
  } |>
  data.table()

real.dt <-
  results |>
  subset(samp.eff == 50) |>
  subset(group.number == 1 & group.rep == 1) |>
  pull(edge.dt) |>
  rbindlist() |>
  subset(type == "real") %>%
  .[,by = .(i,j,n,samp.eff,n.rep,type),
    .(weight = mean(weight))
  ]

SimuNet.dt <-
  results |>
  subset(samp.eff == 50) |>
  subset(group.number == 1 & group.rep == 1) |>
  pull(edge.dt) |>
  rbindlist() |>
  subset(type == "SimuNet") %>%
  .[,by = .(i,j,n,samp.eff,n.rep,type),
    .(weight = mean(weight))
  ]

## plots ----
results.single <-
  results |>
  subset(samp.eff == 50) |>
  subset(group.number == 1 & group.rep == 1) |>
  pull(edge.dt) |>
  rbindlist()
results.single |>
  ggplot(aes(weight, colour = type, fill = type))+
  facet_wrap(. ~ paste0("edge ",i,"-",j),nrow = n.single,scales = "free_y")+
  geom_histogram(binwidth = 1,alpha = 0.15,position = "identity",colour = NA)+
  geom_histogram(data = results.single[type %in% c("SimuNet","real")],
                 binwidth = 1,alpha = 0.15,position = "identity")+
  geom_vline(data = SimuNet.dt,aes(xintercept = weight),
             colour = "tomato",lty = "dashed",size = 1.1)+
  geom_vline(data = real.dt,aes(xintercept = weight),
             colour = "royalblue",lty = "dashed",size = 1.1)+
  geom_point(data = Adj.obs.dt,aes(x = weight,y = 0),fill = "white",
             colour = "tomato",shape = 21,stroke = 1.5,size = 2)+
  scale_colour_manual(values = c("tomato","royalblue",gg_color_hue(6)[2:6]))+
  scale_fill_manual(values = c("tomato","royalblue",gg_color_hue(6)[2:6]),labels = type.str)+
  labs(fill = "Network type",)+
  guides(colour = "none")+
  xlab("Edge weight")+ylab("Count")+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  theme_minimal(16)+theme(legend.position = "top")

results.single |>
  ggplot(aes(weight, colour = type, fill = type))+
  facet_grid(i ~ j,scales = "free_y")+
  geom_density(bw = 1,alpha = 0.25,position = "identity",colour = NA)+
  geom_density(data = results.single[type %in% c("SimuNet","real")],
               bw = 1,alpha = 0.15,position = "identity")+
  geom_vline(data = SimuNet.dt,aes(xintercept = weight),
             colour = "tomato",lty = "dashed",size = 1.1)+
  geom_vline(data = real.dt,aes(xintercept = weight),
             colour = "royalblue",lty = "dashed",size = 1.1)+
  geom_point(data = Adj.obs.dt,aes(x = weight,y = 0),fill = "white",
             colour = "tomato",shape = 21,stroke = 1.5,size = 2)+
  scale_colour_manual(values = c("tomato","royalblue",gg_color_hue(6)[2:6]))+
  scale_fill_manual(values = c("tomato","royalblue",gg_color_hue(6)[2:6]),labels = type.str)+
  labs(fill = "Network type",)+
  guides(colour = "none")+
  xlab("Edge weight")+ylab("Frequency")+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  theme_minimal(16)+theme(legend.position = "top")

results.single |>
  ggplot(aes(weight, colour = type))+
  facet_grid(i ~ j,scales = "free_y")+
  stat_ecdf()+
  stat_ecdf(data = results.single[type %in% c("SimuNet","real")],size = 1.2)+
  geom_vline(data = SimuNet.dt,aes(xintercept = weight),
             colour = "tomato",lty = "dashed",size = 1.1)+
  geom_vline(data = real.dt,aes(xintercept = weight),
             colour = "royalblue",lty = "dashed",size = 1.1)+
  geom_point(data = Adj.obs.dt,aes(x = weight,y = 0),fill = "white",
             colour = "tomato",shape = 21,stroke = 1.5,size = 2)+
  scale_colour_manual(values = c("tomato","royalblue",gg_color_hue(6)[2:6]))+
  labs(colour = "Network type",)+
  xlab("Edge weight")+ylab("Cumulative probability")+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  theme_minimal(16)+theme(legend.position = "top")

results$edge.distance.dt[order(mean.diff)] |>
  dplyr::mutate(type = paste0(y)) |>
  select(c("type","i","j","mean.diff","KS.stat","KS.p","KL","JS","EMD.e","EMD.m")) |>
  GGally::ggpairs(mapping = aes(colour = type,alpha = 0.2),columns = c("mean.diff","KS.stat","KS.p","KL","JS","EMD.e","EMD.m"))+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  theme_minimal()

# Managing batches ----
results |>
  pull(edge.distance.dt) |>
  rbindlist() %>%
  .[,.N,by = .(i,j,samp.eff,y)]

results.multi$edge.distance.dt |>
  # subset(i == 4 & j == 7 ) |>
  ggplot(aes(x = KL,y = y,colour = y,fill = y))+
  facet_grid(i ~ j)+
  geom_density_ridges(aes(group = interaction(i,j,group.rep,y)),jittered_points = TRUE,alpha = 0.03)+
  geom_density_ridges(aes(group = interaction(i,j,group.rep,y)),fill = NA,alpha = 0.8,
                      jittered_points = TRUE,point_shape = "|",point_size = 2.5,position = position_points_jitter(height = 0))+
  # geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_colour_manual(values = paste0(gg_color_hue(5),"50"))+
  scale_y_discrete(labels = x.labs)+
  scale_x_continuous(limits = c(0,NA))+
  ylab("")+xlab("Kullback-Leibler divergence")+
  guides(colour = "none",fill = "none")+
  theme_minimal(12)

results.multi$edge.distance.dt |>
  subset(i <= 3 & j <= 3) |>
  ggplot(aes(x = JS,y = y,colour = as.factor(samp.eff),fill = as.factor(samp.eff)))+
  facet_grid(i ~ j)+
  geom_density_ridges(aes(group = interaction(i,j,samp.eff,y)),alpha = 0.1,scale = .9)+
  geom_density_ridges(aes(group = interaction(i,j,samp.eff,y)),fill = NA,alpha = 0.3,scale = .9,
                      jittered_points = TRUE,point_shape = "|",point_size = 2.5,
                      position = position_points_jitter(height = 0))+
  # geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_colour_manual(values = paste0(gg_color_hue(8),"50"))+
  scale_fill_manual(values = paste0(gg_color_hue(8),"50"))+
  scale_y_discrete(labels = x.labs)+
  scale_x_continuous(limits = c(0,1))+
  ylab("")+xlab("Jensen-Shannon distance")+
  guides(fill = "none")+
  theme_minimal(12)

## Plotting results ----
x.labs <- paste0("... ",c("the real network","another network (similar generation)","an Erdős–Rényi graph (p = 0.5)","a Random network (fixed prob.)","a Random network (variable prob.)"))

KS.stat <-
  results.multi$edge.distance.dt |>
  subset(i <= 3 & j <= 3) |>
  ggplot(aes(x = KS.stat,y = y,colour = y,fill = y))+
  facet_grid(i ~ j)+
  geom_density_ridges(aes(group = interaction(i,j,group.rep,y)),jittered_points = TRUE,alpha = 0.01,scale = .9)+
  geom_density_ridges(aes(group = interaction(i,j,group.rep,y)),fill = NA,alpha = 0.3,scale = .9,
                      jittered_points = TRUE,point_shape = "|",point_size = 2.5,
                      position = position_points_jitter(height = 0))+
  # geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_colour_manual(values = paste0(gg_color_hue(5),"13"))+
  scale_y_discrete(labels = x.labs)+
  scale_x_continuous(limits = c(0,NA))+
  ylab("")+xlab("Kolmogorov-Smirnoff statistic")+
  guides(colour = "none",fill = "none")+
  theme_minimal(12)

KS.p <-
  results.multi$edge.distance.dt |>
  subset(i <= 3 & j <= 3) |>
  ggplot(aes(x = KS.p,y = y,colour = y,fill = y))+
  facet_grid(i ~ j)+
  geom_density_ridges(aes(group = interaction(i,j,group.rep,y)),jittered_points = TRUE,alpha = 0.01,scale = .9)+
  geom_density_ridges(aes(group = interaction(i,j,group.rep,y)),fill = NA,alpha = 0.3,scale = .9,
                      jittered_points = TRUE,point_shape = "|",point_size = 2.5,
                      position = position_points_jitter(height = 0))+
  # geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_colour_manual(values = paste0(gg_color_hue(5),"13"))+
  scale_y_discrete(labels = x.labs)+
  scale_x_continuous(limits = c(0,NA))+
  ylab("")+xlab("Kolmogorov-Smirnoff p-value")+
  guides(colour = "none",fill = "none")+
  theme_minimal(12)

KL <-
  results.multi$edge.distance.dt |>
  subset(i <= 3 & j <= 3) |>
  ggplot(aes(x = KL,y = y,colour = y,fill = y))+
  facet_grid(i ~ j)+
  geom_density_ridges(aes(group = interaction(i,j,group.rep,y)),jittered_points = TRUE,alpha = 0.01,scale = .9)+
  geom_density_ridges(aes(group = interaction(i,j,group.rep,y)),fill = NA,alpha = 0.3,scale = .9,
                      jittered_points = TRUE,point_shape = "|",point_size = 2.5,
                      position = position_points_jitter(height = 0))+
  # geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_colour_manual(values = paste0(gg_color_hue(5),"13"))+
  scale_y_discrete(labels = x.labs)+
  scale_x_continuous(limits = c(0,NA))+
  ylab("")+xlab("Kullback-Leibler divergence")+
  guides(colour = "none",fill = "none")+
  theme_minimal(12)

JS <-
  results.multi$edge.distance.dt |>
  subset(i <= 3 & j <= 3) |>
  ggplot(aes(x = JS,y = y,colour = y,fill = y))+
  facet_grid(i ~ j)+
  geom_density_ridges(aes(group = interaction(i,j,group.rep,y)),jittered_points = TRUE,alpha = 0.01,scale = .9)+
  geom_density_ridges(aes(group = interaction(i,j,group.rep,y)),fill = NA,alpha = 0.3,scale = .9,
                      jittered_points = TRUE,point_shape = "|",point_size = 2.5,
                      position = position_points_jitter(height = 0))+
  # geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_colour_manual(values = paste0(gg_color_hue(5),"13"))+
  scale_y_discrete(labels = x.labs)+
  scale_x_continuous(limits = c(0,NA))+
  ylab("")+xlab("Jensen-Shannon distance")+
  guides(colour = "none",fill = "none")+
  theme_minimal(12)

EMD.e <-
  results.multi$edge.distance.dt |>
  subset(i <= 3 & j <= 3) |>
  ggplot(aes(x = EMD.e,y = y,colour = y,fill = y))+
  facet_grid(i ~ j)+
  geom_density_ridges(aes(group = interaction(i,j,group.rep,y)),jittered_points = TRUE,alpha = 0.01,scale = .9)+
  geom_density_ridges(aes(group = interaction(i,j,group.rep,y)),fill = NA,alpha = 0.3,scale = .9,
                      jittered_points = TRUE,point_shape = "|",point_size = 2.5,
                      position = position_points_jitter(height = 0))+
  # geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_colour_manual(values = paste0(gg_color_hue(5),"13"))+
  scale_y_discrete(labels = x.labs)+
  scale_x_continuous(limits = c(0,NA))+
  ylab("")+xlab("Earth Mover Distance (Euclidean dist.)")+
  guides(colour = "none",fill = "none")+
  theme_minimal(12)

EMD.m <-
  results.multi$edge.distance.dt |>
  subset(i <= 3 & j <= 3) |>
  ggplot(aes(x = EMD.m,y = y,colour = y,fill = y))+
  facet_grid(i ~ j)+
  geom_density_ridges(aes(group = interaction(i,j,group.rep,y)),jittered_points = TRUE,alpha = 0.01,scale = .9)+
  geom_density_ridges(aes(group = interaction(i,j,group.rep,y)),fill = NA,alpha = 0.3,scale = .9,
                      jittered_points = TRUE,point_shape = "|",point_size = 2.5,
                      position = position_points_jitter(height = 0))+
  # geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_colour_manual(values = paste0(gg_color_hue(5),"13"))+
  scale_y_discrete(labels = x.labs)+
  scale_x_continuous(limits = c(0,NA))+
  ylab("")+xlab("Earth Mover Distance (Manhattan dist.)")+
  guides(colour = "none",fill = "none")+
  theme_minimal(12)

gridExtra::grid.arrange(KS.stat,KS.p,KL,JS,EMD.e,EMD.m,ncol = 2,
                        top = grid::textGrob("SimuNet compared to...",
                                             gp=grid::gpar(fontsize = 16,fontface = "bold")))



results.multi$edge.distance.dt %>%
  .[,by = .(y,n,samp.eff,n.rep,group.rep),
     .(
       mean.diff = median(mean.diff),
       mean.diff.min = quantile(mean.diff,probs = 0.025),
       mean.diff.max = quantile(mean.diff,probs = 0.975),
       KS.stat   = median(KS.stat),
       KS.stat.min = quantile(KS.stat,probs = 0.025),
       KS.stat.max = quantile(KS.stat,probs = 0.975),
       KS.p      = median(KS.p),
       KS.p.min = quantile(KS.p,probs = 0.025),
       KS.p.max = quantile(KS.p,probs = 0.975),
       KL        = median(KL),
       KL.min = quantile(KL,probs = 0.025),
       KL.max = quantile(KL,probs = 0.975),
       JS        = median(JS),
       JS.min = quantile(JS,probs = 0.025),
       JS.max = quantile(JS,probs = 0.975),
       EMD.e     = median(EMD.e),
       EMD.e.min = quantile(EMD.e,probs = 0.025),
       EMD.e.max = quantile(EMD.e,probs = 0.975),
       EMD.m     = median(EMD.m),
       EMD.m.min = quantile(EMD.m,probs = 0.025),
       EMD.m.max = quantile(EMD.m,probs = 0.975)
     )
  ] |>
  ggplot(aes(y,KL,colour = y,fill = y))+
  geom_point(alpha = 0.2,position = position_jitterdodge(jitter.width = .5,dodge.width = 2,seed = 42))+
  geom_errorbar(aes(ymin = KL.min, ymax = KL.max),alpha = 0.2,position = position_jitterdodge(jitter.width = .5,dodge.width = 2,seed = 42),width = 0.1)+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  scale_x_discrete(labels = x.labs)+
  xlab("")+ylab("Kullback-Leibler divergence")+
  guides(colour = "none",fill = "none")+
  # coord_flip()+
  theme_minimal(12)




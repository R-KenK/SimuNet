# Load packages and custom functions ----
source(".WIP/experiments_on_modelling_networks_generated_by_complex_processes_tools.R")

# Variables and parameters list ----
set.seed(42)

n.rep   <- 105L
n.group <- 7L
n.each  <- 35L

## generating the list of parameters ----
param.n <-
  # seq(6,30,by = 3)
  c(5,8,10,15,18)
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
  )
param.list <- param.list[sample(1:nrow(param.list))] # perhaps shuffling the row could yield better ETAs

dist.param <-
  param.list |>
  group_by(n,netgen_name,samp.eff,group.number) |>
  summarise() |>
  ungroup() |>
  relocate(c(netgen_name,n,samp.eff)) |>
  data.table() |>
  subset(group.number %in% 1:4)
  # subset(group.number %in% 5:7)
dist.param

# Running simulations ----
## Generating the weighted adjacency matrices ----
start.time <- Sys.time()
start.time
run_simulations(param.list,n.cores = 7)
end.time <- Sys.time()
end.time
end.time - start.time

# Calculating network differences ----
## Preparing to measure edge weights distribution distances ----
prepare_for_distances(dist.param)
Sys.time()

## Measure edge weights distribution distances ----
measure_distances(dist.param)
end.time <- Sys.time()
end.time
end.time - start.time


# Reimporting results ----
query_edgeDistanceDT(edgeDistanceDT.path = ".WIP/simulation.data/edgeDTbis/") |>
  filter(n == 5) |> collect()
arrow::open_dataset(sources = ".WIP/simulation.data/edgeDT/")

## Plotting distances ----
### Single case as boxplots ----
KS.stat <-
  query_edgeDistanceDT() |>
  filter(n == 8 & samp.eff == 75) |>
  filter(group.number == 1) |>
  collect() |>
  plot_distance(dist = "KS.stat",ylab = "Kolmogorov-Smirnoff statistic")

KS.p <-
  query_edgeDistanceDT() |>
  filter(n == 8 & samp.eff == 75) |>
  filter(group.number == 1) |>
  collect() |>
  plot_distance(dist = "KS.p",ylab = "Kolmogorov-Smirnoff p-value")

KL <-
  query_edgeDistanceDT() |>
  filter(n == 8 & samp.eff == 75) |>
  filter(group.number == 1) |>
  collect() |>
  plot_distance(dist = "KL",ylab = "Kullback-Leibler divergence")

JS <-
  query_edgeDistanceDT() |>
  filter(n == 8 & samp.eff == 75) |>
  filter(group.number == 1) |>
  collect() |>
  plot_distance(dist = "JS",ylab = "Jensen-Shannon distance")

# EMD.e <-
#   results |>
#   subset(samp.eff == 75) |>
#   subset(group.number == 1) |>
#   extract_column("edge.distance.dt") |>
#   plot_distance(dist = "EMD.e",ylab = "Earth Mover Distance (Euclidean dist.)")
#
# EMD.m <-
#   results |>
#   subset(samp.eff == 75) |>
#   subset(group.number == 1) |>
#   extract_column("edge.distance.dt") |>
#   plot_distance(dist = "EMD.m",ylab = "Earth Mover Distance (Manhattan dist.)")

gridExtra::grid.arrange(KS.stat,KS.p,KL,JS,#EMD.e,EMD.m,
                        ncol = 2,top = grid::textGrob("SimuNet compared to...",
                                             gp=grid::gpar(fontsize = 16,fontface = "bold")))


### Multiple cases as densities ----
KS.stat <-
  query_edgeDistanceDT() |>
  filter(n == 8 & samp.eff == 75) |>
  collect() |>
  plot_distance(dist = "KS.stat",ylab = "Kolmogorov-Smirnoff statistic",
                geom = "density",x.lims = c(0,1))

KS.p <-
  query_edgeDistanceDT() |>
  filter(n == 8 & samp.eff == 75) |>
  collect() |>
  plot_distance(dist = "KS.p",ylab = "Kolmogorov-Smirnoff p-value",
                geom = "density",x.lims = c(0,1))

KL <-
  query_edgeDistanceDT() |>
  filter(n == 8 & samp.eff == 75) |>
  collect() |>
  plot_distance(dist = "KL",ylab = "Kullback-Leibler divergence",
                geom = "density",x.lims = c(0,NA))

JS <-
  query_edgeDistanceDT() |>
  filter(n == 8 & samp.eff == 75) |>
  collect() |>
  plot_distance(dist = "JS",ylab = "Jensen-Shannon distance",
                geom = "density",x.lims = c(0,1)) # base 2 log make the boundaries [0,1]

# EMD.e <-
#   query_edgeDistanceDT() |>
#   filter(n == 8 & samp.eff == 75) |>
#   collect() |>
#   plot_distance(dist = "EMD.e",ylab = "Earth Mover Distance (Euclidean dist.)",
#                 geom = "density",x.lims = c(0,NA))

# EMD.m <-
#   query_edgeDistanceDT() |>
#   filter(n == 8 & samp.eff == 75) |>
#   collect() |>
#   plot_distance(dist = "EMD.m",ylab = "Earth Mover Distance (Manhattan dist.)",
#                 geom = "density",x.lims = c(0,NA))

gridExtra::grid.arrange(KS.stat,KS.p,KL,JS,# EMD.e,EMD.m,
                        ncol = 2,top = grid::textGrob("SimuNet compared to...",
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


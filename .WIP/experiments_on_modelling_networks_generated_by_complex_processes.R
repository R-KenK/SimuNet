# Load packages and custom functions ----
source(".WIP/experiments_on_modelling_networks_generated_by_complex_processes_tools.R")

# Variables and parameters list ----

n.rep   <- 105L
n.group <- 21L
n.each  <- 56L

## generating the list of parameters ----
param.n <-
  seq(6,24,by = 3)
  # c(5,8,10,15,18)
param.samp.eff <-
  c(c(1,2.5,5) %o% 10^(1:3))[1:7]
  # c(10,25,50,75,100,150,250,500)

param.netgen <-
  list(
    # data.table(netgen_name = "HRCMV",netgen_fun = list(highRankChooseMostValuable_generate_asso)),
    data.table(netgen_name = "GWC",netgen_fun = list(gregWithinClique_generate_asso))
  )

set.seed(42)
param.list <-
  generate_paramList(param.n = param.n,
                     param.samp.eff = param.samp.eff,
                     param.netgen = param.netgen,
                     n.rep = n.rep,
                     n.group = n.group,
                     n.each = n.each
  )[]
# param.list <- param.list[sample(1:nrow(param.list))] # perhaps shuffling the row could yield better ETAs
param.list

# Running simulations ----
## Generating the weighted adjacency matrices ----
start.time <- Sys.time()
start.time
run_simulations(param.list = param.list,delete.tmp = FALSE,n.cores = 7)
end.time <- Sys.time()
end.time
end.time - start.time
# Time difference of 1.649793 hours for param.list of 40 combination of parameters,
# 21 groups with 56 rep each (47040 * 105 simus)

# Calculating network differences ----
## Preparing to measure edge weights distribution distances ----
start.time <- Sys.time()
start.time
prepare_for_distances(param.list = param.list,n.each = n.each,n.chunks = 14)
end.time <- Sys.time()
end.time
end.time - start.time

## Measure edge weights distribution distances ----
start.time <- Sys.time()
start.time
measure_distances(param.list = param.list,n.each = n.each,n.chunks = 7)
end.time <- Sys.time()
end.time
end.time - start.time

# Network metrics bench ----
start.time <- Sys.time()
start.time
calculate_nodeMetrics(param.list = param.list,n.cores = 7)
end.time <- Sys.time()
end.time
end.time - start.time

query_nodeDT() |>
  # filter(group.number == 1) |>
  collect()
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
  filter(n == 9 & samp.eff == 250) |>
  collect() |>
  plot_distance(dist = "KS.stat",ylab = "Kolmogorov-Smirnoff statistic",
                geom = "density",.group = "interaction(group.number,type)",x.lims = c(0,1),.alpha = 0.01)

KS.p <-
  query_edgeDistanceDT() |>
  filter(n == 9 & samp.eff == 250) |>
  collect() |>
  plot_distance(dist = "KS.p",ylab = "Kolmogorov-Smirnoff p-value",
                geom = "density",.group = "interaction(group.number,type)",x.lims = c(0,1),.alpha = 0.01)

KL <-
  query_edgeDistanceDT() |>
  filter(n == 9 & samp.eff == 250) |>
  collect() |>
  plot_distance(dist = "KL",ylab = "Kullback-Leibler divergence",
                geom = "density",.group = "interaction(group.number,type)",x.lims = c(0,NA),.alpha = 0.01)

JS <-
  query_edgeDistanceDT() |>
  filter(n == 9 & samp.eff == 250) |>
  collect() |>
  plot_distance(dist = "JS",ylab = "Jensen-Shannon distance",
                geom = "density",.group = "interaction(group.number,type)",x.lims = c(0,1),.alpha = 0.01) # base 2 log make the boundaries [0,1]

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

p.grid <- gridExtra::grid.arrange(KS.stat,KS.p,KL,JS,# EMD.e,EMD.m,
                        ncol = 2,top = grid::textGrob("How far from the real network are...",
                                             gp=grid::gpar(fontsize = 16,fontface = "bold")))
p.grid
ggsave(filename = ".WIP/simulation.data/Distances densities.png",
       plot = p.grid,width = 15,height = 10,units = "in",dpi = 200)

# Graphic exploration ---------------------------------------------------------------------------------------
## Sampling effort ----
query_edgeDT() |>
  filter(n == 6 & group.number == 1) |>
  collect() |>
  {\(.) .[,weight.scaled := weight / samp.eff]}() |>
  {\(.) .[
    ,
    by = .(netgen_name,n,samp.eff,group.number,type,i,j),

    .(
      weight = median(weight.scaled),
      low = quantile(weight.scaled,probs = 0.025),
      up = quantile(weight.scaled,probs = 0.975)
    )
  ]}() |>
  subset(type %in% c("real","real.bis","SimuNet")) |>
  ggplot(aes(samp.eff,weight,colour = type,fill = type,linetype = type))+
  facet_grid(i ~ j)+
  geom_ribbon(aes(ymin = low,ymax = up),alpha = 0.2)+
  geom_line()+
  geom_point()+
  scale_fill_manual(values = c("tomato","royalblue","darkgreen",rep("grey50",4)))+
  scale_colour_manual(values = c("tomato","royalblue","darkgreen",rep("grey50",4)))+
  theme_bw(16)

## single sim ----
type.str <- c("the real network",
              "the real network bis",
              "SimuNet network",
              "another network (similar generation)",
              "an Erdős–Rényi graph (p = 0.5)",
              "a Random network (fixed prob.)",
              "a Random network (variable prob.)")
n.single <-
  query_edgeDT() |>
  filter(n == 5 & samp.eff == 75) |>
  filter(group.number == 1 & group.rep == 1) |>
  collect() |>
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

